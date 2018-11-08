#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
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
        char* ln = malloc(sizeof(char) * (strlen(line) + 1));
        ln[strlen(line)] = '\0';
        memcpy(ln, line, sizeof(char) * strlen(line));
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
        free(ln);
        free(parts);
      }
    } else { // non-header line
      break;
    }
  }

  return header;
}

bed_line_t *bed_read_line(bed_file_t* bed) {
  //const char* bed_format = "%s\t%d\t%d\t%s\t%d\t%c";

  // bed line parts:
  char chrom[10], name[100];
  int st, en, score;
  char strand;

  bed_line_t *line = malloc(sizeof(bed_line_t));
  
  bed->cur_row++;

  // parse BED line by tokenizing - allows arbitrary fields after BED-6
  char* s = malloc(sizeof(char)*100);
  const char delim[3] = "\t\n";
  size_t len;
  size_t read = getline(&s, &len, bed->fp);
  if(read == -1) return NULL;
  char* token;
  token = strtok(s, delim);
  int i;
  for(i = 0; token != NULL; i++) {
    switch(i) {
      case 0:
        strcpy(line->chrom, token);
        line->chrom[strlen(token)] = NULL; // manually add null-terminator
        break;
      case 1:
        line->st = atoi(token);
        break;
      case 2:
        line->en = atoi(token);
        break;
      case 3:
        strcpy(line->name, token);
        line->name[strlen(token)] = NULL; // manually add null-terminator
        break;
      case 4:
        line->score = atoi(token);
        break;
      case 5:
        line->strand = token[0];
        break;
    }
    token = strtok(NULL, delim);
  }
  /*
  int nfields = fscanf(bed->fp, bed_format, line->chrom, &(line->st), &(line->en), line->name, &(line->score), &(line->strand));
  if(nfields < 6) { // line too short or small or -1 for EOF
    if(nfields == -1) {
      return NULL;
    } else {
      fprintf(stderr, "BED line %d does not appear to be in expected format - STOPPED READING HERE\n", bed->cur_row);
      return NULL;
    }
  }
  */
  free(s);

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

int bed_close(bed_file_t* bed) {
  fclose(bed->fp);
  if(bed->header.sequences != NULL) {
    free(bed->header.sequences);
  }
  if(bed->header.sequence_lengths != NULL) {
    free(bed->header.sequence_lengths);
  }
}
