
#include "colorspace.h"

#define COLOR_0_BASE_ENCODED 'A'
#define COLOR_1_BASE_ENCODED 'T'
#define COLOR_2_BASE_ENCODED 'C'
#define COLOR_3_BASE_ENCODED 'G'

const int MAXLINE = 1024;

int  fasta_to_colorspace(char* fasta_filename, char* output_dir){
	// open the input fasta file
	FILE *input_fasta_file = fopen(fasta_filename, "r");
	if (input_fasta_file == NULL) {
		printf("Error al abrir el fichero %s\n", fasta_filename);	
		return -1;
	}
	
	// create output directory, if necessary
	if (!exists(output_dir)) {
		printf("El directorio '%s' no existe. Creando directorio ...\n", output_dir);
		create_directory(output_dir);
	} else {
		printf("El directorio '%s' ya existe, no es necesario volver a crearlo.\n", output_dir);
	}
	output_dir = strcat(output_dir, "/");

	// fichero de salida
	char *output_file_name = strcat(output_dir, fasta_filename);
	printf("Output file name: %s\n", output_file_name);
	FILE *output_fasta_file = fopen(output_file_name, "w");

	char line[MAXLINE];
	long num_lineas = 0;
	char *old_line;
	char last_nucleotide_of_line;
	char colorspace_line[MAXLINE];
	while (fgets(line, MAXLINE, input_fasta_file)) {
		if (line[0] == '>') {
			fprintf(output_fasta_file, "%s", line);
			last_nucleotide_of_line = NULL;
		} else {
			last_nucleotide_of_line = nucleotide_sequence_to_color_space(last_nucleotide_of_line, line, colorspace_line);
			fprintf(output_fasta_file, "%s\n", colorspace_line);
			//fputs(colorspace_line, output_fasta_file);
			//free(colorspace_line);
		}
		num_lineas++;
	}
	fclose(input_fasta_file);
	fclose(output_fasta_file);

	printf("Numero de lineas leidas: %ld\n", num_lineas);

	return 0;
}

char nucleotide_sequence_to_color_space(char last_nucleotide_of_line, char *nucleotide_line, char *colorspace_line) {
	char color;
	int contador = 0;
	//colorspace_line = malloc(MAXLINE);
	if (last_nucleotide_of_line != NULL) {
		//printf("Primera linea del fasta\n");
		colorspace_line[contador] = dinucleotide_to_color( last_nucleotide_of_line, nucleotide_line[0]);
		contador ++;
	}
	int seq_length = strlen(nucleotide_line) - 1;
	int i;
	for (i=0; i<seq_length-1; i++) {
		colorspace_line[contador] = dinucleotide_to_color( nucleotide_line[i], nucleotide_line[i+1]);
		contador++;
	}
	colorspace_line[contador] = 0;
	last_nucleotide_of_line = nucleotide_line[seq_length-1];
	return last_nucleotide_of_line;
}

char dinucleotide_to_color(char n1, char n2) {
	// TODO: seguro que hay una forma mas rapida y menos cutre
	char color;
	if (n1 == 'N') { n1 = 'A'; }
	if (n2 == 'N') { n2 = 'A'; }
	if (n1 == 'A' && n2 == 'A') {
		color = COLOR_0_BASE_ENCODED;
	} else if (n1 == 'A' && n2 == 'C') {
		color = COLOR_1_BASE_ENCODED;
	} else if (n1 == 'A' && n2 == 'G') {
		color = COLOR_2_BASE_ENCODED;
	} else if (n1 == 'A' && n2 == 'T') {
		color = COLOR_3_BASE_ENCODED;
	} else if (n1 == 'C' && n2 == 'A') {
		color = COLOR_1_BASE_ENCODED;
	} else if (n1 == 'C' && n2 == 'C') {
		color = COLOR_0_BASE_ENCODED;
	} else if (n1 == 'C' && n2 == 'G') {
		color = COLOR_3_BASE_ENCODED;
	} else if (n1 == 'C' && n2 == 'T') {
		color = COLOR_2_BASE_ENCODED;
	} else if (n1 == 'G' && n2 == 'A') {
		color = COLOR_2_BASE_ENCODED;
	} else if (n1 == 'G' && n2 == 'C') {
		color = COLOR_3_BASE_ENCODED;
	} else if (n1 == 'G' && n2 == 'G') {
		color = COLOR_0_BASE_ENCODED;
	} else if (n1 == 'G' && n2 == 'T') {
		color = COLOR_1_BASE_ENCODED;
	} else if (n1 == 'T' && n2 == 'A') {
		color = COLOR_3_BASE_ENCODED;
	} else if (n1 == 'T' && n2 == 'C') {
		color = COLOR_2_BASE_ENCODED;
	} else if (n1 == 'T' && n2 == 'G') {
		color = COLOR_1_BASE_ENCODED;
	} else if (n1 == 'T' && n2 == 'T') {
		color = COLOR_0_BASE_ENCODED;
	}
	return color;
}

int  cs_fastq_to_base_encoding(char* input_fastq_filename, char* output_fastq_filename){
	// open the input and output fastq files
	printf("Creating input readers ...\n");
	fastq_file_t *input_fastq_file = fastq_fopen_mode(input_fastq_filename, "r");
	fastq_file_t *output_fastq_file = fastq_fopen_mode(output_fastq_filename, "w");

	// read every read from the input fastq file
	printf("Reading sequences from input fastq file ...\n");
	fastq_read_t *input_read = fastq_read_new("","",""),
				 *output_read;
	char bs_encoded_sequence[MAXLINE];
	while (fastq_fread(input_read, input_fastq_file)) {
		cs_sequence_to_base_space_encoding(input_read->sequence, bs_encoded_sequence);
		output_read = fastq_read_new(input_read->id, bs_encoded_sequence, input_read->quality);
		fastq_fwrite(output_read, 1, output_fastq_file);
		// TODO: borrar las reads, ¿es necesario?
		//fastq_read_free(output_read);
	}

	// close the files
	printf("Closing fastq files ...\n");
	fastq_fclose(input_fastq_file);
	fastq_fclose(output_fastq_file);
}

void cs_sequence_to_base_space_encoding(char* cs_sequence, char* bs_encoded_cs_sequence) {
	// transform the encoding of each char in the input sequence
	int seq_length = strlen(cs_sequence);
	int i;
	for (i=0; i<seq_length; i++) {
		bs_encoded_cs_sequence[i] = base_space_color_encoding(cs_sequence[i]);
	}
	bs_encoded_cs_sequence[i] = 0;
}

char base_space_color_encoding(char color) {
	char base;
	if (color == '0') {
		base = COLOR_0_BASE_ENCODED;
	} else if (color == '1') {
		base = COLOR_1_BASE_ENCODED;
	} else if (color == '2') {
		base = COLOR_2_BASE_ENCODED;
	} else if (color == '3') {
		base = COLOR_3_BASE_ENCODED;
	}
	return base;
}
