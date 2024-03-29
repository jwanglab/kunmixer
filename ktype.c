#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "incl/klib/khash.h"
#include "incl/klib/ksort.h"
#include "kseq_by_type.h"
#include "incl/klib/kvec.h"
#include "bed.h"
#include <zlib.h>
#include <bzlib.h>


// have to reorder params to make this work with kseq
int fileread(FILE* f, char* buffer, int size) {
  return fread(buffer, 1, size, f);
}

// creates string:[array of uint8] hash
// to map chrom names to sequences
KHASH_MAP_INIT_STR(refSeq, char*);

// creates string:uint32 hash
// to map chrom names to lengths
KHASH_MAP_INIT_STR(refLen, uint32_t);

// creates string:int hash
// to map chrom names to IDs (tid)
KHASH_MAP_INIT_STR(refId, int);


// creates int:uint32 hash
// for mapping loci within a chromosome to an INDEX into the locus vector
KHASH_MAP_INIT_INT(lociMap, uint32_t);

// creates string:lociMap hash
// for mapping chromosome to the corresponding locus map
KHASH_MAP_INIT_STR(chromMap, khash_t(lociMap)*);

// A,C,G,T,N,Unknown
typedef uint16_t* locus;

typedef struct locus_allele {
  uint32_t locus_id;
  uint8_t allele_id;
  int8_t offset; // position of the variable allele in the kmer key
} locus_allele;

// creates string:locus_allele hash
// for mapping kmers to the corresponding locus id and allele
KHASH_MAP_INIT_STR(kmerMap, locus_allele);

const char ALLELES[] = {'A', 'C', 'G', 'T', 'N'};
const int N_ALLELES = 5;


void add_with_error(khash_t(kmerMap) *kmer_hash, char* kmer, int k, int locus_id, int snp_offset) {
  int i, absent;
  khint_t bin;
  char* dup;
  for(i = 0; i < k; i++) {
    if(i == (-1*snp_offset)) continue;
    char tmp = kmer[i];
    int alt_id = 0;
    for(alt_id = 0; alt_id < N_ALLELES; alt_id++) {
      if(tmp == ALLELES[alt_id]) continue;
      kmer[i] = ALLELES[alt_id];
      dup = malloc(sizeof(char) * (strlen(kmer) + 1));
      dup[strlen(kmer)] = '\0';
      memcpy(dup, kmer, sizeof(char) * strlen(kmer));
      bin = kh_put(kmerMap, kmer_hash, dup, &absent);
      if(!absent) {
        fprintf(stderr, "kmer '%s' already exists - maybe from a different locus?!\n", kmer);
        continue;
      }
      locus_allele lal = {locus_id, alt_id, snp_offset};
      /*
      int j;
      for(j = 0; j < k+snp_offset; j++) printf(" ");
      printf("%s     with SNP at %d and err at %d is locus %d allele %d\n", kmer, snp_offset, i, lal.locus_id, lal.allele_id);
      */
      kh_val(kmer_hash, bin) = lal;
    }
    kmer[i] = tmp;
  }
}


char compl(char a) {
  switch(a) {
    case 'A':
      return 'T';
      break;
    case 'C':
      return 'G';
      break;
    case 'G':
      return 'C';
      break;
    case 'T':
      return 'A';
      break;
  }
  return 'N';
}

int bz2read_wrapper(BZFILE *b, void *buf, int len) {
  int bzerror;
  int l = BZ2_bzRead ( &bzerror, b, buf, len );
  return l;
}

// init kseq struct
KSEQ_INIT_READER(gzFile, gzread, gz)
KSEQ_INIT_READER(BZFILE*, bz2read_wrapper, bz2)

int main(int argc, char *argv[]) {

  if(argc < 4) {
    fprintf(stderr, "Usage: ktype <k> <BED> <reference FASTA> <read FASTA/Q> [<read FASTA/Q> ...]\n");
    fprintf(stderr, "Not enough arguments.\n");
    return 1;
  }
  int k = atoi(argv[1]);
  char *bed_file = argv[2];
  char *ref_fasta = argv[3];
  int n_read_files = argc - 4;
  char **read_files = malloc(sizeof(char*) * n_read_files);
  int i;
  for(i = 4; i < argc; i++) {
    read_files[i-4] = argv[i];
  }

  khint_t bin, subbin; // hash bin (result of kh_put)
  int absent;

  // ordered vector of loci in the BED file
  kvec_t(locus) alleles;
  kv_init(alleles);
  kvec_t(bed_line_t*) loci;
  kv_init(loci);
  kvec_t(char) loci_ref_allele;
  kv_init(loci_ref_allele);


  // load ref FASTA file
  //
  khash_t(refSeq) *ref = kh_init(refSeq);
  khash_t(refLen) *rlen = kh_init(refLen);
  khash_t(refId) *rid = kh_init(refId);

  gzFile gzfp;
  kseq_t_gz *seq;
  int l;

  gzfp = gzopen(ref_fasta, "r");
  if(!gzfp) {
    fprintf(stderr, "File '%s' not found\n", ref_fasta);
    return 1;
  }
  fprintf(stderr, "Reading fasta file: %s\n", ref_fasta);
  seq = kseq_init_gz(gzfp);

  int refid = 0;
  char* dup;
  while ((l = kseq_read_gz(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    //printf("Reading %s (%i bp).\n", seq->name.s, l);

    dup = malloc(sizeof(char) * (strlen(seq->name.s) + 1));
    dup[strlen(seq->name.s)] = '\0';
    memcpy(dup, seq->name.s, sizeof(char) * strlen(seq->name.s));

    // seq array
    bin = kh_put(refSeq, ref, dup, &absent);
    // copy the seq read from kseq to a new heap here - this is pretty fast and the easiest way to implement right now (see kseq.h)
    kh_val(ref, bin) = malloc(sizeof(char)*l);
    memcpy(kh_val(ref, bin), seq->seq.s, sizeof(char)*l);

    // sequence length
    bin = kh_put(refLen, rlen, dup, &absent);
    kh_val(rlen, bin) = l;

    // ref ID (indexed order in FASTA)
    bin = kh_put(refId, rid, dup, &absent);
    kh_val(rid, bin) = refid++;
  }

  gzclose(gzfp);
  kseq_destroy_gz(seq);


  // load BED file
  //
  fprintf(stderr, "Reading BED file '%s'\n", bed_file);
  bed_file_t bed = bed_init(bed_file);
  if(bed.fp == NULL) {
    return 1;
  }
  khash_t(chromMap) *regions = kh_init(chromMap); // our nested hash map (l0: chrom, l1: pos)
  bed_line_t *entry = bed_read_line(&bed);

  khash_t(kmerMap) *kmer_hash = kh_init(kmerMap);
  char* kmer = malloc(sizeof(char)*(k+1)); // this will be reused repeatedly
  kmer[k] = '\0';
  char* kmer_rc = malloc(sizeof(char)*(k+1));
  kmer_rc[k] = '\0';

  int tot_kmer = 0;
  int dup_kmer = 0;

  while(entry != NULL) {
    //fprintf(stderr, "%s %d %d\n", entry->chrom, entry->st, entry->en);
    kv_push(bed_line_t*, loci, entry);
    
    //fprintf(stderr, "chrom lookup\n");
    // look up chromosome in first-level hash
    dup = malloc(sizeof(char) * (strlen(entry->chrom) + 1));
    dup[strlen(entry->chrom)] = '\0';
    memcpy(dup, entry->chrom, sizeof(char) * strlen(entry->chrom));
    bin = kh_put(chromMap, regions, dup, &absent);
    if(absent) {
      kh_val(regions, bin) = kh_init(lociMap);
    }
    //fprintf(stderr, "second lookup\n");
    // look up locus in second-level hash
    subbin = kh_put(lociMap, kh_val(regions, bin), entry->st, &absent);
    int locus_id = kv_size(alleles);
    kh_val(kh_val(regions, bin), subbin) = locus_id;
    locus l = calloc(6, sizeof(uint16_t));
    kv_push(locus, alleles, l);

    //fprintf(stderr, "ref lookup\n");
    bin = kh_get(refId, rid, entry->chrom);
    absent = (bin == kh_end(rid)); 
    if(absent) {
      fprintf(stderr, "Chromosome '%s' not found in reference\n", entry->chrom);
      return 1;
    }
    int tid = kh_val(rid, bin);
    //fprintf(stderr, "tid: %d\n", tid);

    // construct k-mers in this region with different alleles
    // and add to k-mer hash
    bin = kh_get(refSeq, ref, entry->chrom);
    char* rs = kh_val(ref, bin);
    // add ref allele to the list
    kv_push(char, loci_ref_allele, rs[entry->st]);
    int offset; // from the target locus position
    // this will crash if a locus is within (k)nt of the start or end of a chromosome
    for(offset = -1*(k-1); offset <= 0; offset++) {
      //fprintf(stderr, "offset: %d\n", offset);
      memcpy(kmer, rs+entry->st+offset, sizeof(char)*k);
      //fprintf(stderr, "kmer: %s\n", kmer);
      //make uppercase
      for(i = 0; i < k; i++) {
        kmer[i] = toupper(kmer[i]);
        kmer_rc[k-1-i] = compl(kmer[i]);
      }
      int alt_id = 0;
      char* dup;
      for(alt_id = 0; alt_id < N_ALLELES; alt_id++) {
        kmer[-1*offset] = ALLELES[alt_id];
        //fprintf(stderr, " mod kmer: %s\n", kmer);
        dup = malloc(sizeof(char) * (strlen(kmer)+1));
        dup[strlen(kmer)] = '\0';
        memcpy(dup, kmer, sizeof(char) * strlen(kmer));
        //fprintf(stderr, " dup: %s\n", dup);
        //fprintf(stderr, " done dup\n");
        bin = kh_get(kmerMap, kmer_hash, dup);
        absent = (bin == kh_end(kmer_hash));
        //fprintf(stderr, " tot_kmer++\n");
        tot_kmer++;
        //fprintf(stderr, " if(!absent)\n");
        if(!absent) {
          //fprintf(stderr, "kmer '%s' already exists - maybe from a different locus?!\n", kmer);
        //fprintf(stderr, " dup_kmer++\n");
          dup_kmer++;
        //fprintf(stderr, " continue\n");
          continue;
        } else {
          bin = kh_put(kmerMap, kmer_hash, dup, &absent);
        }
        //fprintf(stderr, " create lal\n");
        locus_allele lal = {locus_id, alt_id, offset};

        //for(i = 0; i < k+offset; i++) fprintf(stderr, " ");
        //fprintf(stderr, "%s     is %d allele %d\n", kmer, lal.locus_id, lal.allele_id);

        //fprintf(stderr, " kh_val(kmer_hash, bin) = lal\n");
        kh_val(kmer_hash, bin) = lal;

        // and reverse-complement
        kmer_rc[k-1+offset] = compl(ALLELES[alt_id]);
        //fprintf(stderr, " rc kmer: %s\n", kmer_rc);
        dup = malloc(sizeof(char) * (strlen(kmer_rc)+1));
        dup[strlen(kmer_rc)] = '\0';
        memcpy(dup, kmer_rc, sizeof(char) * strlen(kmer_rc));
        //fprintf(stderr, " dup: %s\n", dup);
        //fprintf(stderr, " done dup\n");
        bin = kh_get(kmerMap, kmer_hash, dup);
        absent = (bin == kh_end(kmer_hash));
        //fprintf(stderr, " tot_kmer++\n");
        tot_kmer++;
        //fprintf(stderr, " if(!absent)\n");
        if(!absent) {
          //fprintf(stderr, "kmer '%s' (-) already exists - maybe from a different locus?!\n", kmer_rc);
        //fprintf(stderr, " dup_kmer++\n");
          dup_kmer++;
        //fprintf(stderr, " continue\n");
          continue;
        } else {
          bin = kh_put(kmerMap, kmer_hash, dup, &absent);
        }
        //fprintf(stderr, " kh_val(kmer_hash, bin) = lal\n");
        kh_val(kmer_hash, bin) = lal;

        // then add all k-mers w/one nt error (except our variant):
        //add_with_error(kmer_hash, kmer, k, locus_id, offset);
        // this adds ~5-10% more hits, but many are wrong, so it's not clear that it's helping
      }
    }

    entry = bed_read_line(&bed);
    //if(kv_size(loci) >= 7) break;
  }
  //bed_close(&bed);
  fprintf(stderr, "%d of %d (%f%%) k-mers duplicated\n", dup_kmer, tot_kmer, (float)dup_kmer/tot_kmer*100);


  void* seq2;
  // go through each read fasta/q file in turn
  int r, s, n; // r: index into read_files, s: index into read sequence, n: number of reads processed
  int hits = 0;
  for(r = 0; r < n_read_files; r++) {
    FILE* f;
    BZFILE* b = NULL;
    int bzerror;

    // init kseq struct
    if(strncmp(read_files[r]+strlen(read_files[r])-3, ".gz", 3) == 0) {
      fprintf(stderr, "File '%s' is gzipped\n", read_files[r]);
      gzfp = gzopen(read_files[r], "r");
      if(!gzfp) {
        fprintf(stderr, "File '%s' not found\n", read_files[r]);
        return 1;
      }
      seq2 = kseq_init_gz(gzfp);
    } else if (strncmp(read_files[r]+strlen(read_files[r])-4, ".bz2", 4) == 0) {
      
      fprintf(stderr, "File '%s' is bzip2'd\n", read_files[r]);
      f = fopen(read_files[r], "r");
      if (!f) {
        fprintf(stderr, "File '%s' not found\n", read_files[r]);
        return 1;
      }
      b = BZ2_bzReadOpen(&bzerror, f, 0, 0, (void *)NULL, 0);
      if(!b) {
        fprintf(stderr, "File '%s' identified as BZ2 but failed to open (?)\n", read_files[r]);
        return 1;
      }
      seq2 = kseq_init_bz2(b);
    }

    //printf("Reading fasta file: %s\n", read_files[r]);

    while ((l = b ? kseq_read_bz2((kseq_t_bz2*)seq2) : kseq_read_gz((kseq_t_gz*)seq2)) > 0) {
      // name: seq->name.s, seq: seq->seq.s, length: l
      //fprintf(stderr, "Reading %s (%i bp): %s (%d) %s (%d).\n", seq2->name.s, l, seq2->seq.s, seq2->seq.l, seq2->qual.s, seq2->qual.l);

      for(s = 0; s < l; s=s+k) {
        if(b)
          memcpy(kmer, ((kseq_t_bz2*)seq2)->seq.s+(s+k >= l ? l-k : s), k);
        else
          memcpy(kmer, ((kseq_t_gz*)seq2)->seq.s+(s+k >= l ? l-k : s), k);
        bin = kh_get(kmerMap, kmer_hash, kmer);
        absent = (bin == kh_end(kmer_hash)); 

        if(!absent) {
          hits++;
          locus_allele lal = kh_val(kmer_hash, bin);
          kv_A(alleles, lal.locus_id)[lal.allele_id]++;
          //for(i = 0; i < k+lal.offset; i++) printf(" ");
          //fprintf(stderr, "%s     at read %d pos %d, hit locus %d allele %d\n", kmer, n, s-lal.offset, lal.locus_id, lal.allele_id);
        }
      }

      n++;
      if(n % 1000 == 0) {
        //break;
      }
    }
    //fprintf(stderr, "Ended reading %s (%i bp): %s.\n", seq2->name.s, l, seq2->seq.s);

    if(b) {
      BZ2_bzReadClose(&bzerror, b);
      fclose(f);
      kseq_destroy_bz2((kseq_t_bz2*)seq2);
    } else {
      gzclose(gzfp);
      kseq_destroy_gz((kseq_t_gz*)seq2);
    }
  }
  //fprintf(stderr, "%d k-mers hit a variant\n", hits);


  // output "variant" positions
  //FILE* fout = fopen()
  printf("chrom\tstart\tend\tref_allele\tA\tC\tG\tT\tN\tUnknown\n");
  for (i = 0; i < kv_size(loci); i++) {
    // at this point, we're assuming these regions are in the hash(es) (and NOT checking)
    bin = kh_get(chromMap, regions, kv_A(loci, i)->chrom);
    subbin = kh_get(lociMap, kh_val(regions, bin), kv_A(loci, i)->st);
    locus al = kv_A(alleles, kh_val(kh_val(regions, bin), subbin));
    printf("%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", kv_A(loci, i)->chrom, kv_A(loci, i)->st, kv_A(loci, i)->en, kv_A(loci_ref_allele, i), al[0], al[1], al[2], al[3], al[4], al[5]);
  }

  // free BED lines and containing vector
  for(i = 0; i < kv_size(loci); i++) {
    free(kv_A(loci, i));
  }
  kv_destroy(loci);

  for(i = 0; i < kv_size(alleles); i++) {
    free(kv_A(alleles, i));
  }
  kv_destroy(alleles);

  free(kmer);

  return 0;
}

