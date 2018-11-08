#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "incl/htslib/htslib/sam.h"
#include "incl/klib/khash.h"
#include "incl/klib/ksort.h"
#include "incl/klib/kseq.h"

/*
 * fastdiff.c
 *
 * Jeremy Wang
 * 20180419
 *
 * As fast as possible compute the loci where the average read from
 * the bam file differs from the reference
*/

// have to reorder params to make this work with kseq
int fileread(FILE* f, char* buffer, int size) {
  return fread(buffer, 1, size, f);
}

// init kseq struct
KSEQ_INIT(FILE*, fileread);

// creates string:[array of uint8] hash
KHASH_MAP_INIT_STR(refSeq, char*);

// creates string:uint32 hash
KHASH_MAP_INIT_STR(refLen, uint32_t);


int main(int argc, char *argv[]) {

  if(argc < 3) {
    printf("Usage: fastdiff <BAM> <reference FASTA>\n");
    printf("Not enough arguments.\n");
    return -1;
  }
  char *bam_file = argv[1];
  char *ref_fasta = argv[2];


  // load FASTA file

  khash_t(refSeq) *ref = kh_init(refSeq);
  khash_t(refLen) *rlen = kh_init(refLen);

  FILE* fp;
  kseq_t *seq, *nextseq;
  int l, absent;
  char* dup;

  fp = fopen(ref_fasta, "r");
  seq = kseq_init(fp);
  printf("Reading fasta file: %s\n", ref_fasta);

  khint_t bin; // hash bin (result of kh_put)

  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    printf("Reading %s (%i bp).\n", seq->name.s, l);

    // seq array
    dup = malloc(sizeof(char) * (strlen(seq->name.s) + 1));
    dup[strlen(seq->name.s)] = '\0';
    memcpy(dup, seq->name.s, sizeof(char) * strlen(seq->name.s));
    bin = kh_put(refSeq, ref, dup, &absent);
    // copy the seq read from kseq to a new heap here - this is pretty fast and the easiest way to implement right now (see kseq.h)
    kh_val(ref, bin) = malloc(sizeof(char)*l);
    memcpy(kh_val(ref, bin), seq->seq.s, sizeof(char)*l);

    // sequence length
    bin = kh_put(refLen, rlen, dup, &absent);
    kh_val(rlen, bin) = l;
  }

  fclose(fp);
  kseq_destroy(seq);


  // load BAM file
  //
  samFile *bam;
  bam_hdr_t *header;
  bam1_t *aln;
  int ret_val;

  bam = sam_open(bam_file, "rb");
  if (bam == NULL) {
    fprintf(stderr, "Error opening \"%s\"\n", bam_file);
    return -1;
  }
  header = sam_hdr_read(bam);
  if (header == NULL) {
    fprintf(stderr, "Couldn't read header for \"%s\"\n", bam_file);
    return -1;
  }
  // construct array from reference information so that we can look it up with read.tid
  char **ref_array = malloc(sizeof(char*) * header->n_targets);
  int8_t **diffs = malloc(sizeof(int8_t*) * header->n_targets); // (ref_id x ref_len) array of variant counts
  int *rlen_array = malloc(sizeof(uint32_t) * header->n_targets); // has to be a plain int because that's what kseq gives out
  //printf("BAM targets:\n");
  int i, j;
  for (i = 0; i < header->n_targets; i++) {
    diffs[i] = calloc(header->target_len[i], sizeof(int8_t));
    //printf("%s (%u bp)\n", header->target_name[i], header->target_len[i]);
    bin = kh_get(refSeq, ref, header->target_name[i]);
    ref_array[i] = kh_value(ref, bin);
    /*
    printf("Seq: ");
    int j;
    for(j = 0; j < 100; j++) {
      printf("%c", ref_array[i][j]);
    }
    printf("\n");
    */
    bin = kh_get(refLen, rlen, header->target_name[i]);
    int fa_len = kh_value(rlen, bin);
    if (fa_len != header->target_len[i]) { // target_len is a uint32_t
      printf("WARNING: Reference fasta length (%i) and BAM length (%u) of %s do not agree.\n", fa_len, header->target_len[i], header->target_name[i]);
    }
    rlen_array[i] = fa_len;
  }

  aln = bam_init1();

  int32_t read_len = 0;

  uint32_t ct = 0;
  while ((ret_val = sam_read1(bam, header, aln)) >= 0) {

    if (aln->core.flag & 4) { // unmapped
      continue;
    }

    int32_t tid = aln->core.tid;
    int32_t qlen = aln->core.l_qseq;
    int32_t pos = aln->core.pos;
    int32_t endpos = bam_endpos(aln) - 1;
    int32_t* cigar = bam_get_cigar(aln); // lower 4 bits are cigar operation, upper 28 are the length
    uint8_t* qseq = bam_get_seq(aln); // 4 bits each (1: A, 2: C, 4: G, 8: T, 15: N)
    int32_t qpos = 0;
    int32_t o;
    uint8_t rev = bam_is_rev(aln);
    //printf("query of len %d aligned to target %d (%s) at pos %d (%c), cigar: ", qlen, tid, header->target_name[tid], pos, (rev ? '-' : '+'));

    for(o = 0; o < aln->core.n_cigar; o++) {
      char opc = bam_cigar_opchr(cigar[o]);
      int oplen = bam_cigar_oplen(cigar[o]);
      //printf("%c%d, ", opc, oplen);
      int32_t qstartpos = qpos;
      int32_t startpos = pos;
      if(opc == 'M') {
        /*
        while(qpos < qstartpos + oplen) {
          char qc = ".AC.G...T......N"[bam_seqi(qseq, qpos)];
          //printf("%c/%c ", qc, ref_array[tid][pos]);
          printf("%d%c ", bam_seqi(qseq, qpos), qc);
          //pos++;
          qpos++;
        }
        printf("\n");
        while(pos < startpos + oplen) {
          printf("%c", ref_array[tid][pos]);
          pos++;
        }
        printf("\n");
        */
        while(qpos < qstartpos + oplen) {
          char qc = ".AC.G...T......N"[bam_seqi(qseq, qpos)];
          char rc = ref_array[tid][pos];
          // convert to uppercase
          if(qc >= 97) qc = qc - 32;
          if(rc >= 97) rc = rc - 32;
          // sometimes the query sequence is messed up - I need to figure that out
          if((qc == 'A' || qc == 'C' || qc == 'G' || qc == 'T') && (rc == 'A' || rc == 'C' || rc == 'G' || rc == 'T')) {
            // set non-over(under)flowing count of alt vs. ref alleles
            if(qc != rc && diffs[tid][pos] < 127) {
              diffs[tid][pos]++;
            } else if(diffs[tid][pos] > -128) {
              diffs[tid][pos]--;
            }
          }
          pos++;
          qpos++;
        }
      } else if(opc == 'I' || opc == 'S') { // insertion or softclip both consume the query only
        qpos = qpos + oplen;
      } else if(opc == 'D') {
        pos = pos + oplen;
      } else if(opc == 'H') { // hard clip, does nothing...
        continue;
      } else {
        printf("WARNING: Unhandled CIGAR op character '%c' at target %d, position %d\n", opc, tid, pos);
        break;
      }
    }

    ct++;
    if(ct % 10000 == 0) {
      printf("%d reads processed\n", ct);
    }

    //printf("\n");
  }

  // output "variant" positions
  //FILE* fout = fopen()
  for (i = 0; i < header->n_targets; i++) {
    for (j = 0; j < header->target_len[i]; j++) {
      if(diffs[i][j] > 1) {
        //fprintf(fout, "%s\t%d\t%d", header->target_name[i], j, diffs[i][j]);
        printf("%s\t%d\t%d\n", header->target_name[i], j, diffs[i][j]);
      }
    }
  }

  // clean up this first pass through the BAM file
  bam_hdr_destroy(header);

  ret_val = sam_close(bam);
  if (ret_val < 0) {
    fprintf(stderr, "Error closing input BAM.\n");
    return -1;
  }

  bam_destroy1(aln);

  return 0;
}

