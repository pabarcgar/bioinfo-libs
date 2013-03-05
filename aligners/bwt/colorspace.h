#ifndef COLORSPACE_H
#define COLORSPACE_H

#include <stdio.h>
#include <string.h>

#include "commons/file_utils.h"

int fasta_to_colorspace(char* fasta_filename, char* output_directory);
char nucleotide_sequence_to_color_space(char last_nucleotide_of_line, char *nucleotide_line, char *colorspace_line);
char dinucleotide_to_color(char n1, char n2);

#endif // COLORSPACE_H
