#include <stdio.h>

#include "colorspace.h"

#define MIN_OPTIONS 3

int main(int argc, char** argv) {
  if (argc < MIN_OPTIONS) {
    printf("Usage: %s Dna_reference_Path Output_dir\n", argv[0]);
    //exit(-1);
    return -1;
  }
  char *fasta_file_name = argv[1];
  char *output_dir = argv[2];

  fasta_to_colorspace(fasta_file_name, output_dir);
}

