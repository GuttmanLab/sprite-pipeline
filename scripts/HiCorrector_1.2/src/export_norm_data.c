/*
 * export_norm_data.c
 *
 *  Created on: March 24, 2015
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
	fprintf(stream, "\texport_norm_data <input raw matrix file> <#rows/columns> <has header line in input file?> <has header column input file?> <memory size (MB)> <input bias vector file> <fixed row sum after normalization> <output normalized matrix file>\n");
//	fprintf(stream, "\nOptions:\n\t<output file format>: FULL_MATRIX_TEXT\n");
	fprintf(stream, "\nVersion: 1.2\n");
	fprintf(stream, "\nAuthor: Wenyuan Li\n");
	fprintf(stream, "Date: March 24, 2015\n");
}

int main(int argc, char *argv[]) {
	int ret, i;
	char input_mat_file[MAXCHAR], input_bias_vec_file[MAXCHAR];
	float mem_size_of_task;
	int part_size;
	float x;
	float row_sum_=1.0;
	float row_sum_after_norm=1.0;
	int total_column;
	int has_header_line;
	int has_header_column;
	char output_file[MAXCHAR];
	VEC_FLOAT bias_vec;
	VEC_INT part_sizes;
	int start_row, end_row;
	MATRIX_FLOAT mat; // submatrix loaded, its size is "4*part_size*total_column" bytes (float is 4 bytes in both 32bit and 64bit machine: see http://docs.oracle.com/cd/E18752_01/html/817-6223/chp-typeopexpr-2.html)

	if (argc!=9) {
		fprintf(stderr,"Error: export_norm_data arguments are insufficient.\n\n");
		print_usage(stderr);
		exit(EXIT_FAILURE);
	}
	strcpy(input_mat_file, argv[1]);
	total_column = atoi(argv[2]);
	has_header_line = atoi(argv[3]);
	has_header_column = atoi(argv[4]);
	mem_size_of_task = (float)atof(argv[5]);
	strcpy(input_bias_vec_file, argv[6]);
	row_sum_after_norm = (float)atof(argv[7]);
	strcpy(output_file, argv[8]);

	printf("export_norm_data arguments:\n");
	printf("\tinput raw matrix file: %s\n", input_mat_file);
	printf("\t\thas header line: %d\n", has_header_line);
	printf("\t\thas header column: %d\n", has_header_column);
	printf("\ttotal_rows: %d\n", total_column);
	printf("\tinput bias vector file: %s\n", input_bias_vec_file);
	printf("\tmemory size for this job: %g MB\n", mem_size_of_task);
	printf("\trow sum after normalization: %f\n", row_sum_after_norm);
	printf("\toutput file: %s\n", output_file);
	fflush( stdout );

	ret = quick_check_partitions_by_mem_size(total_column, total_column, mem_size_of_task, &x);
	if (ret==FALSE) {
		fprintf(stderr, "Error: max_memory_size (%gMB) for this task cannot afford to load the matrix (size=%d rows X %d columns, each row as a partition needs %lf MB memory)\n", mem_size_of_task, total_column, total_column, x);
		exit(EXIT_FAILURE);
	}

	calc_partitions_by_mem_size(total_column, total_column, mem_size_of_task, &part_sizes);
	part_size = max_VEC_INT(part_sizes);
	printf("All data is divided to %d parts: ", part_sizes.n);
	put_VEC_INT(stdout, part_sizes, ", ");
	printf("\npartition sizes (%d partitions): ", part_sizes.n);
	put_VEC_INT(stdout, part_sizes, ", ");
	printf("\n");

	init_VEC_FLOAT( &bias_vec );
	ret = read_VEC_FLOAT( &bias_vec, input_bias_vec_file );
	if (ret==0) {
		fprintf(stderr, "Error: Wrong reading %s\n", input_bias_vec_file);
		exit(EXIT_FAILURE);
	}
	init_MATRIX_FLOAT( &mat );

	for (i=0, start_row=1; i<part_sizes.n; i++) {
		// the i-th part
		// start_row and end_row are 1-based index
		end_row = start_row + part_sizes.v[i] - 1;
		printf("\tPart %d:\t%d\t%d\n", i+1, start_row, end_row);
		ret = read_efficient_MATRIX_FLOAT_with_header(input_mat_file, start_row, end_row, total_column, has_header_line, has_header_column, &mat);
		if (ret==EXIT_SUCCESS) {
			// update matrix: B = diag(1./d_save)*A*diag(1./d_save);
			elementwise_div_vector_MATRIX_FLOAT( &mat, bias_vec);
			if (i==0) {
				row_sum_ = norm1_float(mat.v[0], mat.total_column) / row_sum_after_norm;
			}
			elementwise_div_value_MATRIX_FLOAT( &mat, row_sum_ );
			// append this matrix to file
			if (i==0)
				write_MATRIX_FLOAT(output_file, mat);
			else
				append_MATRIX_FLOAT_text(output_file, mat);
		} else {
			fprintf(stderr, "Error: Wrong reading on lines %d - %d\n",
					start_row, end_row);
			exit(EXIT_FAILURE);
		}
		start_row = end_row + 1;
	}
	// free space
	free_VEC_INT( &part_sizes );
	free_MATRIX_FLOAT( &mat );
	free_VEC_FLOAT( &bias_vec );
	printf("Done.\n");
	return EXIT_SUCCESS;
}

