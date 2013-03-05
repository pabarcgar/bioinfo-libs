#include <stdio.h>

#include "colorspace.h"

#define MIN_OPTIONS 3

int main(int argc, char** argv) {
  if (argc < MIN_OPTIONS) {
    printf("Usage: %s input_fastq_filename output_fastq_filename\n", argv[0]);
    //exit(-1);
    return -1;
  }
  char *input_fastq_file_name = argv[1];
  char *output_fastq_file_name = argv[2];

  cs_fastq_to_base_encoding(input_fastq_file_name, output_fastq_file_name);
}
