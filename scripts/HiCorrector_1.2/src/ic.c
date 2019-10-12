/*
 * ic.c
 *
 *  Created on: Aug 4, 2014
 *  Revised on: Oct 17, 2014
 *		1. Add the function that can read input matrix file with header line or column

 *      Author: Wenyuan Li
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "matutils.h"

void print_usage(FILE* stream)
{
	fprintf(stream, "Usage:\n");
	fprintf(stream, "\tic <input matrix file> <#rows/columns> <max iterations> <has header line in input file?> <has header column in input file?> <output file>\n");
	fprintf(stream, "\n\t\thas header line: '1' indicates input matrix file has a header line; otherwise 0;");
	fprintf(stream, "\n\t\thas header column: '1' indicates input matrix file has a header column; otherwise 0.\n");
	fprintf(stream, "\nVersion: 1.1\n");
	fprintf(stream, "Author: Wenyuan Li\n");
	fprintf(stream, "Date: October 17, 2014\n");
}

int main(int argc, char *argv[]) {
	int ret;
	time_t t_start, t_end;
	double t_elapse; // seconds
	char input_mat_file[MAXCHAR];
	float p = 1.0;
	int total_column;
	int maxiter;
	int has_header_line;
	int has_header_column;
	char output_file[MAXCHAR];
	MATRIX_FLOAT mat; // data matrix loaded
	int iter;
	VEC_FLOAT b, t;

	if (argc!=7) {
		fprintf(stderr,"Error: Arguments are insufficient.\n\n");
		print_usage(stderr);
		exit(EXIT_FAILURE);
	}
	strcpy(input_mat_file, argv[1]);
	total_column = atoi(argv[2]);
	maxiter = atoi(argv[3]);
	has_header_line = atoi(argv[4]);
	has_header_column = atoi(argv[5]);
	strcpy(output_file, argv[6]);

	printf("ic arguments:\n");
	printf("\tinput file: %s\n", input_mat_file);
	printf("\t\thas header line: %d\n", has_header_line);
	printf("\t\thas header column: %d\n", has_header_column);
	printf("\ttotal_column: %d\n", total_column);
	printf("\tmaxiter: %d\n", maxiter);
	printf("\toutput file: %s\n", output_file);
	fflush( stdout );

	t_start = time(NULL);
	create_ones_VEC_FLOAT( total_column, &b);
	create_VEC_FLOAT( total_column, &t);
	init_MATRIX_FLOAT( &mat );
	ret = read_efficient_MATRIX_FLOAT_with_header(input_mat_file, 1, total_column, total_column, has_header_line, has_header_column, &mat);
	printf("Load data done\n");
	fflush( stdout );
	if (ret==EXIT_SUCCESS) {
		for (iter=0; iter<maxiter; iter++) {
			printf("\tIteration %d:\n", iter+1); fflush( stdout );
			// get new norms from updated matrix: t = normp(O,p);
			norm_rows_MATRIX_FLOAT(mat, p, &t);
			// update matrix: O = diag(1./t)*O*diag(1./t);
			elementwise_div_vector_MATRIX_FLOAT(&mat, t);
			// update the bias vector: b = b.*t;
			elementwise_multi_vector_VEC_FLOAT(&b, 1, total_column, t);

			t_end = time(NULL);
			t_elapse = difftime( t_end, t_start );
			printf("\t\ttime - %d:%d:%d\n", (int)floor(t_elapse/3600.0), (int)floor(fmod(t_elapse,3600.0)/60.0), (int)fmod(t_elapse,60.0));
			fflush( stdout );
		}
	} else {
		fprintf(stderr, "Error: Wrong reading file '%s'\nExit.\n", input_mat_file);
		exit(EXIT_FAILURE);
	}

	t_end = time(NULL);
	t_elapse = difftime( t_end, t_start );
	printf(">> ic computation total time - %d:%d:%d\n", (int)floor(t_elapse/3600.0), (int)floor(fmod(t_elapse,3600.0)/60.0), (int)fmod(t_elapse,60.0));
	fflush( stdout );

	// write output
	t_start = time(NULL);
	printf("Write output to file: %s\n", output_file);
	write_VEC_FLOAT(b, output_file);
	t_end = time(NULL);
	t_elapse = difftime( t_end, t_start );
	printf("Write time: %d:%d:%d\n", (int)floor(t_elapse/3600.0), (int)floor(fmod(t_elapse,3600.0)/60.0), (int)fmod(t_elapse,60.0));

	// free space
	free_MATRIX_FLOAT( &mat );
	free_VEC_FLOAT( &b );
	free_VEC_FLOAT( &t );
	printf("Done.\n");
	return EXIT_SUCCESS;
}
