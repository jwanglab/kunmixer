
typedef struct bed_header_t {
  char* version;
  size_t num_sequences;
  char** sequences;
  size_t* sequence_lengths;
} bed_header_t;

typedef struct bed_file_t {
  FILE *fp;
  bed_header_t header;
  size_t cur_row;
} bed_file_t;

typedef struct bed_line_t {
  char chrom[10];
  int st;
  int en;
  char name[100];
  int score;
  char strand;
} bed_line_t;

bed_header_t parse_header(FILE* bed_fp);
bed_line_t *bed_read_line(bed_file_t* bed);
bed_file_t bed_init(char* f);
