#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "incl/htslib/htslib/sam.h"
#include "incl/htslib/htslib/khash.h"
#include "incl/htslib/htslib/ksort.h"
#include "incl/htslib/htslib/kseq.h"
#include "incl/klib/kvec.h"
#include "bed.h"

// wikipedia
void show_bits(unsigned int x) {
    int i; 
    for(i=(sizeof(int)*8)-1; i>=0; i--)
      (x&(1u<<i))?putchar('1'):putchar('0');
    printf("\n");
}

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
// to map chrom names to sequences
KHASH_MAP_INIT_STR(refSeq, char*);

// creates string:uint32 hash
// to map chrom names to lengths
KHASH_MAP_INIT_STR(refLen, uint32_t);

// creates string:int hash
// to map chrom names to IDs (tid)
KHASH_MAP_INIT_STR(refId, int);


// creates int:uint8 hash
// for mapping loci within a chromosome to an INDEX into the locus vector
KHASH_MAP_INIT_INT(lociMap, uint8_t);

// creates string:lociMap hash
// for mapping chromosome to the corresponding locus map
KHASH_MAP_INIT_STR(chromMap, khash_t(lociMap)*);

typedef struct locus{
  uint16_t a;
  uint16_t c;
  uint16_t g;
  uint16_t t;
  uint16_t n;
  uint16_t u; // unknown (sometimes '.')
} locus;


int main(int argc, char *argv[]) {

  if(argc < 4) {
    fprintf(stderr, "Usage: fasttype <BAM> <BED> <reference FASTA>\n");
    fprintf(stderr, "Not enough arguments.\n");
    return -1;
  }
  char *bam_file = argv[1];
  char *bed_file = argv[2];
  char *ref_fasta = argv[3];

  khint_t bin, subbin; // hash bin (result of kh_put)
  int absent;

  // ordered vector of loci in the BED file
  kvec_t(locus) alleles;
  kv_init(alleles);
  kvec_t(bed_line_t*) loci;
  kv_init(loci);
  kvec_t(char) loci_ref_allele;
  kv_init(loci_ref_allele);


  // load FASTA file
  //
  khash_t(refSeq) *ref = kh_init(refSeq);
  khash_t(refLen) *rlen = kh_init(refLen);
  khash_t(refId) *rid = kh_init(refId);

  FILE* fp;
  kseq_t *seq, *nextseq;
  int l;

  fp = fopen(ref_fasta, "r");
  seq = kseq_init(fp);
  //printf("Reading fasta file: %s\n", ref_fasta);

  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    //printf("Reading %s (%i bp).\n", seq->name.s, l);

    // seq array
    bin = kh_put(refSeq, ref, strdup(seq->name.s), &absent);
    // copy the seq read from kseq to a new heap here - this is pretty fast and the easiest way to implement right now (see kseq.h)
    kh_val(ref, bin) = malloc(sizeof(char)*l);
    memcpy(kh_val(ref, bin), seq->seq.s, sizeof(char)*l);

    // sequence length
    bin = kh_put(refLen, rlen, strdup(seq->name.s), &absent);
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
  int *rlen_array = malloc(sizeof(uint32_t) * header->n_targets); // has to be a plain int because that's what kseq gives out
  //printf("BAM targets (%d):\n", header->n_targets);
  int i, j;
  for (i = 0; i < header->n_targets; i++) {
    //printf("%s (%u bp)\n", header->target_name[i], header->target_len[i]);
    bin = kh_get(refSeq, ref, header->target_name[i]);
    ref_array[i] = kh_value(ref, bin);
    bin = kh_get(refLen, rlen, header->target_name[i]);
    int fa_len = kh_value(rlen, bin);
    if (fa_len != header->target_len[i]) { // target_len is a uint32_t
      fprintf(stderr, "WARNING: Reference fasta length (%i) and BAM length (%u) of %s do not agree.\n", fa_len, header->target_len[i], header->target_name[i]);
    }
    rlen_array[i] = fa_len;
    bin = kh_put(refId, rid, header->target_name[i], &absent);
    kh_val(rid, bin) = i;
  }


  // initialize empty locus mask to mark loci in the BED file
  // 
  uint8_t **locus_mask = malloc(sizeof(uint8_t*) * header->n_targets);
  for(i = 0; i < header->n_targets; i++) {
    locus_mask[i] = calloc((header->target_len[i]/8) + 1, sizeof(uint8_t)); // one bit per locus
  }


  // load BED file
  //
  bed_file_t bed = bed_init(bed_file);
  khash_t(chromMap) *regions = kh_init(chromMap); // our nested hash map (l0: chrom, l1: pos)
  bed_line_t *entry = bed_read_line(&bed);
  while(entry != NULL) {
    //printf("%s %d %d\n", entry->chrom, entry->st, entry->en);
    kv_push(bed_line_t, loci, entry);
    
    // look up chromosome in first-level hash
    bin = kh_put(chromMap, regions, strdup(entry->chrom), &absent);
    if(absent) {
      kh_val(regions, bin) = kh_init(lociMap);
    }
    // look up locus in second-level hash
    subbin = kh_put(lociMap, kh_val(regions, bin), entry->st, &absent);
    kh_val(kh_val(regions, bin), subbin) = kv_size(alleles);
    locus l = {0,0,0,0,0,0};
    kv_push(locus, alleles, l);

    // set locus mask
    bin = kh_get(refId, rid, entry->chrom);
    // TODO: check to make sure bin exists - this segfaults if the BAM file was broken somehow
    int tid = kh_val(rid, bin);
    locus_mask[tid][entry->st/8] |= (1<<(entry->st%8));

    // add ref allele to the list
    kv_push(char, loci_ref_allele, ref_array[tid][entry->st]);

    entry = bed_read_line(&bed);
  }


  // process reads...
  //
  aln = bam_init1();

  int32_t read_len = 0;

  uint32_t ct = 0;
  while ((ret_val = sam_read1(bam, header, aln)) >= 0) {

    ct++;
    if(ct % 1000000 == 0) {
      //printf("%d reads processed\n", ct);
    }

    if (aln->core.flag & 4) { // unmapped
      continue;
    }

    int32_t tid = aln->core.tid;
    int32_t qlen = aln->core.l_qseq;
    int32_t pos = aln->core.pos;
    int32_t endpos = bam_endpos(aln) - 1;

    //printf("trying %d %d %d\n", tid, pos, endpos);

    if(tid >= header->n_targets) {
      continue;
    }

    // check if any position covered by this read is in our loci set
    // the fast way - using the locus mask
    uint8_t overlaps = 0;
    for(i = pos/8; i < endpos/8 + 1; i++) {
      if(locus_mask[tid][i] > 0) {
        overlaps = 1;
        break;
      }
    }
    if(!overlaps) {
      continue;
    }

    //printf("%d %d-%d overlaps\n", tid, pos, endpos);

    int32_t* cigar = bam_get_cigar(aln); // lower 4 bits are cigar operation, upper 28 are the length
    uint8_t* qseq = bam_get_seq(aln); // 4 bits each (1: A, 2: C, 4: G, 8: T, 15: N)
    int32_t qpos = 0;
    int32_t o;
    uint8_t rev = bam_is_rev(aln);
    //printf("query of len %d aligned to target %d (%s) at pos %d (%c), cigar: ", qlen, tid, header->target_name[tid], pos, (rev ? '-' : '+'));
 
    // look up nested hashes (slower than mask lookup)
    bin = kh_get(chromMap, regions, strdup(header->target_name[tid]));
    absent = (bin == kh_end(regions)); 
    if(absent) {
      continue;
    }
    khash_t(lociMap) *lmap = kh_val(regions, bin);

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
          //printf("%d%c ", bam_seqi(qseq, qpos), qc);
          printf("%c", qc);
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
            // look up locus in second-level hash
            subbin = kh_get(lociMap, lmap, pos);
            absent = (subbin == kh_end(lmap)); 
            if(!absent) {
              if(qc == 'A')
                kv_A(alleles, kh_val(lmap, subbin)).a++;
              else if(qc == 'C')
                kv_A(alleles, kh_val(lmap, subbin)).c++;
              else if(qc == 'G')
                kv_A(alleles, kh_val(lmap, subbin)).g++;
              else if(qc == 'T')
                kv_A(alleles, kh_val(lmap, subbin)).t++;
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
        fprintf(stderr, "WARNING: Unhandled CIGAR op character '%c' at target %d, position %d\n", opc, tid, pos);
        break;
      }
    }
  }

  // output "variant" positions
  //FILE* fout = fopen()
  printf("chrom\tstart\tend\tref_allele\tA\tC\tG\tT\tN\tUnknown\n");
  for (i = 0; i < kv_size(loci); i++) {
    // at this point, we're assuming these regions are in the hash(es) (and NOT checking)
    bin = kh_get(chromMap, regions, strdup(kv_A(loci, i)->chrom));
    subbin = kh_get(lociMap, kh_val(regions, bin), kv_A(loci, i)->st);
    locus al = kv_A(alleles, kh_val(kh_val(regions, bin), subbin));
    printf("%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", kv_A(loci, i)->chrom, kv_A(loci, i)->st, kv_A(loci, i)->en, kv_A(loci_ref_allele, i), al.a, al.c, al.g, al.t, al.n, al.u);
  }

  // clean up locus mask
  for(i = 0; i < header->n_targets; i++) {
    free(locus_mask[i]);
  }
  free(locus_mask);

  // clean up this first pass through the BAM file
  bam_hdr_destroy(header);

  ret_val = sam_close(bam);
  if (ret_val < 0) {
    fprintf(stderr, "Error closing input BAM.\n");
    return -1;
  }

  bam_destroy1(aln);

  // free BED lines and containing vector
  for(i = 0; i < kv_size(loci); i++) {
    free(kv_A(loci, i));
  }
  kv_destroy(loci);

  kv_destroy(alleles);

  return 0;
}

