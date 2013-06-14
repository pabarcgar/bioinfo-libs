#include <stdio.h>

#include "colorspace.h"

#define MIN_OPTIONS 3

void fastq_nt_to_ntcs(char *in_filename, char *out_filename);

int main(int argc, char** argv) {
  if (argc < MIN_OPTIONS) {
    printf("Usage: %s Dna_reference_Path Output_dir\n", argv[0]);
    //exit(-1);
    return -1;
  }
  char *fasta_file_name = argv[1];
  char *output_dir = argv[2];

  //  fasta_to_colorspace(fasta_file_name, output_dir);
  fastq_nt_to_ntcs(argv[1], argv[2]);
}

void fastq_nt_to_ntcs(char *in_filename, char *out_filename) {

  printf("Creating input readers ...\n");
  fastq_file_t *input_fastq_file = fastq_fopen_mode(in_filename, "r");
  fastq_file_t *output_fastq_file = fastq_fopen_mode(out_filename, "w");
  
  // read every read from the input fastq file
  printf("Reading sequences from input fastq file ...\n");
  fastq_read_t fq_read;
  char *cs_string;
  while (fastq_fread(&fq_read, input_fastq_file)) {
    cs_string = (char *) calloc(strlen(fq_read.sequence), sizeof(char));

    nucleotide_sequence_to_color_space(NULL, fq_read.sequence, cs_string);
    fq_read.quality[strlen(cs_string)] = 0;

    if (fq_read.sequence != NULL) free(fq_read.sequence);

    fq_read.sequence = cs_string;
    fastq_fwrite(&fq_read, 1, output_fastq_file);
    
    // free memory
    if (fq_read.id != NULL) free(fq_read.id);
    if (fq_read.sequence != NULL) free(fq_read.sequence);
    if (fq_read.quality != NULL) free(fq_read.quality);
  }
  
  // close the files
  printf("Closing fastq files ...\n");
  fastq_fclose(input_fastq_file);
  fastq_fclose(output_fastq_file);
}

