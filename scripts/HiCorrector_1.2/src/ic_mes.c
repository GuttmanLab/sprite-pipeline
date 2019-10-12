/*
 ============================================================================
 Name        : ic_mes.c
   Created on: Jun 25, 2014
   Revised on: Oct 17, 2014
     1. Add the function that can read input matrix file with header line or column

 Author      : Wenyuan Li
 Version     : 1.1

 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "matutils.h"

void print_usage(FILE* stream)
{
	fprintf(stream, "Usage:\n");
	fprintf(stream, "\tic_mes <input matrix file> <memory size (MB)> <#rows/columns> <max iterations> <has header line in input file?> <has header column in input file?> <output file>\n");
	fprintf(stream, "\nVersion: 1.1\n");
	fprintf(stream, "Author: Wenyuan Li\n");
	fprintf(stream, "Date: October 17, 2014\n");
}

int main(int argc, char *argv[]) {
	int ret, i;
	time_t t_start, t_end;
	double t_elapse; // seconds
	char input_mat_file[MAXCHAR];
	float p = 1.0;
	float mem_size_of_task;
	int part_size;
	int total_column;
	int maxiter;
	int has_header_line;
	int has_header_column;
	char output_file[MAXCHAR];
	VEC_INT part_sizes;
	int start_row, end_row, iter;
	MATRIX_FLOAT mat; // submatrix loaded, its size is "4*part_size*total_column" bytes (float is 4 bytes in both 32bit and 64bit machine: see http://docs.oracle.com/cd/E18752_01/html/817-6223/chp-typeopexpr-2.html)
	MATRIX_FLOAT d_save; // accumulated factors (2 vectors, one for used, the other for backup), its size is"4*2*total_column" bytes
	VEC_FLOAT d_temp; // norms of rows, its size is "4*part_size" bytes
	float *d_save_use_ptr, *d_save_new_ptr;
	float x;

	if (argc!=8) {
		fprintf(stderr,"Error: ic_mes arguments are insufficient.\n\n");
		print_usage(stderr);
		exit(EXIT_FAILURE);
	}
	strcpy(input_mat_file, argv[1]);
	mem_size_of_task = (float)atof(argv[2]);
	total_column = atoi(argv[3]);
	maxiter = atoi(argv[4]);
	has_header_line = atoi(argv[5]);
	has_header_column = atoi(argv[6]);
	strcpy(output_file, argv[7]);

	printf("ic_mes arguments:\n");
	printf("\tinput file: %s\n", input_mat_file);
	printf("\t\thas header line: %d\n", has_header_line);
	printf("\t\thas header column: %d\n", has_header_column);
	printf("\tmemory size for this job: %g MB\n", mem_size_of_task);
	printf("\ttotal_rows: %d\n", total_column);
	printf("\tmaxiter: %d\n", maxiter);
	printf("\toutput file: %s\n", output_file);
	fflush( stdout );

	ret = quick_check_partitions_by_mem_size(total_column, total_column, mem_size_of_task, &x);
	if (ret==FALSE) {
		fprintf(stderr, "Error: max_memory_size (%gMB) for this task cannot afford to compute on the matrix (size=%d rows X %d columns, each row as a partition needs %lf MB memory)\n", mem_size_of_task, total_column, total_column, x);
		exit(EXIT_FAILURE);
	}

	t_start = time(NULL);
	calc_partitions_by_mem_size(total_column, total_column, mem_size_of_task, &part_sizes);
	part_size = max_VEC_INT(part_sizes);
	printf("All data is divided to %d parts: ", part_sizes.n);
	put_VEC_INT(stdout, part_sizes, ", ");
	printf("\npartition sizes (%d partitions): ", part_sizes.n);
	put_VEC_INT(stdout, part_sizes, ", ");
	printf("\n");

	init_MATRIX_FLOAT( &mat );
	create_ones_VEC_FLOAT( part_size, &d_temp );
	create_ones_MATRIX_FLOAT( 2, total_column, &d_save);

	for (iter=0; iter<maxiter; iter++) {
		// when iter is even, use d_save[0] to compute and deposit the updated to d_save[1]
		// when iter is odd,  use d_save[1] to compute and deposit the updated to d_save[0]
		if (iter%2==0) {
			d_save_use_ptr = d_save.v[0];
			d_save_new_ptr = d_save.v[1];
		} else {
			d_save_use_ptr = d_save.v[1];
			d_save_new_ptr = d_save.v[0];
		}
		printf("\tIteration %d:\n", iter+1);
		for (i=0, start_row=1; i<part_sizes.n; i++) {
			// the i-th part
			// start_row and end_row are 1-based index
			end_row = start_row + part_sizes.v[i] - 1;
			printf("\t\tPart %d:\t%d\t%d\n", i+1, start_row, end_row);
			ret = read_efficient_MATRIX_FLOAT_with_header(input_mat_file, start_row, end_row, total_column, has_header_line, has_header_column, &mat);
			if (ret==EXIT_SUCCESS) {
				// update matrix: B = diag(1./d_save)*A*diag(1./d_save);
				if (iter!=0)
					elementwise_div_vector_MATRIX_FLOAT_2( &mat, d_save_use_ptr, total_column);
				// get new norms from updated matrix: d_temp = normp(B,p);
				norm_rows_MATRIX_FLOAT(mat, p, &d_temp);
				// update accumulated factors: d_save = d_save.*d_temp;
				elementwise_multi_vector_VEC_FLOAT_2(d_save_use_ptr, start_row, end_row, d_temp, d_save_new_ptr);
			} else {
				fprintf(stderr, "Error: Wrong reading on lines %d - %d\n",
						start_row, end_row);
				exit(EXIT_FAILURE);
			}
			start_row = end_row + 1;
		}
	}

	t_end = time(NULL);
	t_elapse = difftime( t_end, t_start );
	printf(">> ic_mes computation total time: %d:%d:%d\n", (int)floor(t_elapse/3600.0), (int)floor(fmod(t_elapse,3600.0)/60.0), (int)fmod(t_elapse,60.0));

	// write output
	t_start = time(NULL);
	printf("Write output to file: %s\n", output_file);
	write_VEC_FLOAT_2(d_save_new_ptr, total_column, "\n", output_file);
	t_end = time(NULL);
	t_elapse = difftime( t_end, t_start );
	printf("Write time: %d:%d:%d\n", (int)floor(t_elapse/3600.0), (int)floor(fmod(t_elapse,3600.0)/60.0), (int)fmod(t_elapse,60.0));

	// free space
	free_VEC_INT( &part_sizes );
	free_MATRIX_FLOAT( &mat );
	free_MATRIX_FLOAT( &d_save );
	free_VEC_FLOAT( &d_temp );
	printf("Done.\n");
	return EXIT_SUCCESS;
}
