#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "bed.h"

bed_header_t parse_header(FILE* bed_fp) {
  bed_header_t header;
  // 1024 maximum ref sequences
  header.num_sequences = 0;
  header.sequences = malloc(sizeof(char*) * 1024);
  header.sequence_lengths = malloc(sizeof(size_t) * 1024);

  char line[1024]; // maximum line size is 1024 chars
  const char* delim = "\t";
  const char* field_delim = ":";
  char *key, *val;
  char *hd_type;
  char **parts;
  
  // get lines as long as they start with '@' (header)
  int c;
  int i, j;
  while(1) {
    c = fgetc(bed_fp); // get first character
    ungetc(c, bed_fp); // QUICK! put it back! (basically just a peek)
    if((char)c == '@') { // header line
      if(fgets(line, sizeof line, bed_fp) != NULL) {
        char* ln = strdup(line);
        parts = malloc(sizeof(char*) * 10); // maximum 10 fields in a single header line
        i = 0;
        parts[i] = strtok(ln, delim);
        while(parts[i] != NULL) {
          i++;
          parts[i] = strtok(NULL, delim);
        }
        int nfields = i;
        for(i = 0; i < nfields; i++) {
          if(strlen(parts[i]) == 0 || (int)parts[i][0] == 13) { // 13: carriage return
            continue;
          }
          if(strcmp(parts[i], "@HD") == 0) {
            hd_type = parts[i];
          } else if(strcmp(parts[i], "@SQ") == 0) {
            hd_type = parts[i];
            header.num_sequences++;
          } else {
            key = strtok(parts[i], field_delim);
            val = strtok(NULL, field_delim);
            if(strcmp(key, "VN") == 0) {
              header.version = val;
            } else if(strcmp(key, "SN") == 0) {
              header.sequences[header.num_sequences-1] = val;
            } else if(strcmp(key, "LN") == 0) {
              header.sequence_lengths[header.num_sequences-1] = atoi(val);
            }
          }
        }
        //for(j = 0; j < nfields; j++) {
        //  free(parts[j]);
        //}
        free(parts);
      }
    } else { // non-header line
      break;
    }
  }

  return header;
}

bed_line_t *bed_read_line(bed_file_t* bed) {
  const char* bed_format = "%s\t%d\t%d\t%s\t%d\t%c";

  // bed line parts:
  char chrom[10], name[100];
  int st, en, score;
  char strand;

  bed_line_t *line = malloc(sizeof(bed_line_t));

  int nfields = fscanf(bed->fp, bed_format, line->chrom, &(line->st), &(line->en), line->name, &(line->score), &(line->strand));
  if(nfields != 6) { // line too short or small or -1 for EOF
    if(nfields == -1) {
      return NULL;
    } else {
      fprintf(stderr, "A line does not appear to be in expected BED format - STOPPED READING HERE\n");
      return NULL;
    }
  }

  return line;
}

bed_file_t bed_init(char* f) {

  // open bed file from path
  FILE *bed_fp = fopen(f, "r");
  if( bed_fp == NULL ) {
    fprintf(stderr, "Error reading bed file '%s'\n", f);
  }

  bed_header_t header;
  /*
  printf("Parsing header\n");
  header = parse_header(bed_fp);

  printf("BED version %s\n", header.version);
  printf("%d sequences found\n", header.num_sequences);
  */

  bed_file_t bed = {bed_fp, header, 0};

  return bed;
}
