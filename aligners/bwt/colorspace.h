#ifndef COLORSPACE_H
#define COLORSPACE_H

#include <stdio.h>
#include <string.h>

#include "commons/file_utils.h"
#include "bioformats/fastq/fastq_file.h"
#include "bioformats/fastq/fastq_read.h"

int fasta_to_colorspace(char* fasta_filename, char* output_filename);
char nucleotide_sequence_to_color_space(char last_nucleotide_of_line, char *nucleotide_line, char *colorspace_line);
char dinucleotide_to_color(char n1, char n2);
int  cs_fastq_to_base_encoding(char* input_fastq_filename, char* output_fastq_filename);
void cs_sequence_to_base_space_encoding(char* cs_sequence);
char base_space_color_encoding(char color);

#endif // COLORSPACE_H
