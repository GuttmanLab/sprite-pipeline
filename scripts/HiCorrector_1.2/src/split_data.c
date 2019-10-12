/*
 * split_data.c
 *
 *  Created on: August 14, 2014
 *      Author: Wenyuan Li
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "matutils.h"

void split_data_and_write_for_one_worker(char *input_mat_file, int taskid, int start_row, int end_row,
	int total_column, float mem_size_of_task, char *output_dir, char *job_id)
{
	int i, ret;
	int cnt_row, sub_start_row, sub_end_row;
	VEC_INT part_sizes; // its size is small and thus ignored.
	char temp_data_file[MAXCHAR];
	// start_row and end_row are 1-based
	// receive info of which partition assigned to this worker task
	cnt_row = end_row - start_row + 1;
	calc_partitions_by_mem_size(cnt_row, total_column, mem_size_of_task, &part_sizes);
	printf("worker %d: %d - %d (%d parts), part sizes:", taskid,start_row,end_row,part_sizes.n);
	put_VEC_INT(stdout, part_sizes, ", ");

	for (i=0, sub_start_row=start_row; i<part_sizes.n; i++) {
		// the i-th part
		// start_row and end_row are 1-based index
		sub_end_row = sub_start_row + part_sizes.v[i] - 1;
		// In iteration 0, extract the partition and write it to the temporary file
		sprintf(temp_data_file,"%s.temp.%s.task%d.subpart%d",output_dir,job_id,taskid,i);
		printf("\tworker %d: - read whole data file '%s', write sub-part to temp file '%s', for sub-part %d (%d, %d)\n", taskid,
				input_mat_file, temp_data_file, i, sub_start_row, sub_end_row); fflush( stdout );
		ret = read_write_part_of_MATRIX_FLOAT(input_mat_file, sub_start_row, sub_end_row, total_column, temp_data_file);
		sub_start_row = sub_end_row + 1;
	} // end of for "iterate on parts"

	// free space
	free_VEC_INT( &part_sizes );
	printf ("Done (worker %d).\n", taskid); fflush( stdout );
}

void print_usage(FILE* stream)
{
	fprintf(stream, "Usage:\n");
	fprintf(stream, "\tsplit_data <input matrix file> <#rows/columns> <output directory> <#task> <memory size (MB) for each task> <job_id>\n");
	fprintf(stream, "\nVersion: 1.0\n");
	fprintf(stream, "\nNote: #task includes both master and worker tasks.\n");
	fprintf(stream, "\nAuthor: Wenyuan Li\n");
	fprintf(stream, "Date: August 14, 2014\n");
}

int main(int argc, char ** argv) {
	int	numtasks,              /* number of tasks in partition */
		taskid,                /* a task identifier */
		numworkers;            /* number of worker tasks */
	int i, ret;
	float x;
	time_t t_start, t_end;
	double t_elapse; // seconds
	char input_mat_file[MAXCHAR];
	int total_column = atoi(argv[2]);
	float mem_size_of_task = (float)atof(argv[5]);
	char output_file[MAXCHAR], output_dir[MAXCHAR], job_id[MAXCHAR];
	VEC_INT part_sizes;
	int start_row, end_row;
	VEC_INT2 part_info_of_workers;

	if (argc!=7) {
		fprintf(stderr,"Error: split_data arguments (%d) are insufficient.\n\n", argc);
		print_usage(stderr);
		exit(EXIT_FAILURE);
	}
	numtasks = atoi(argv[4]);
	strcpy(input_mat_file, argv[1]);
	strcpy(output_dir, argv[3]);
	strcpy(job_id, argv[6]);
	i = strlen(output_dir);
	if (output_dir[i-1]!='/')
		strcat(output_dir,"/");

    numworkers = numtasks-1;

    /**************************** master task ************************************/
	printf("split_data arguments:\n");
	printf("\tinput file: %s\n", input_mat_file);
	printf("\ttotal_column: %d\n", total_column);
	printf("\toutput file: %s\n", output_file);
	printf("\t#task: %d\n", numtasks);
	printf("\tmemory size for each cpu: %g\n", mem_size_of_task);
	fflush( stdout );

	t_start = time(NULL);

	partitions2( total_column, numworkers, &part_sizes ); // #rows for each task
	printf("split_data has %d tasks (both master and %d workers).\n", numtasks, numworkers);
	printf("master: partition sizes (%d partitions): ", part_sizes.n);
	put_VEC_INT(stdout, part_sizes, ", ");
	fflush( stdout );
	if (numworkers!=part_sizes.n) {
		fprintf(stderr, "Error (master): #parts %d != numworkers %d", part_sizes.n, numworkers);
		exit(EXIT_FAILURE);
	}

	// send partition info to each worker
	printf("split data\n"); fflush( stdout );
	create_VEC_INT2(numworkers, &part_info_of_workers);
	for (taskid=1, start_row=1; taskid<=numworkers; taskid++) {
		// the i-th part
		// start_row and end_row are 1-based index
		i = taskid - 1;
		end_row = start_row + part_sizes.v[i] - 1;
		ret = quick_check_partitions_by_mem_size(part_sizes.v[i], total_column, mem_size_of_task, &x);
		if (ret==FALSE) {
			fprintf(stderr, "Error: max_memory_size (%gMB) for each task cannot afford to compute on the matrix (size=%d rows X %d columns, each row as a partition needs %lf MB memory)\n", mem_size_of_task, part_sizes.v[i], total_column, x);
			exit(EXIT_FAILURE);
		}
		part_info_of_workers.v1[taskid-1] = start_row;
		part_info_of_workers.v2[taskid-1] = end_row;
		printf("master: worker %d - (%d, %d)\n", taskid, start_row, end_row); fflush( stdout );

		split_data_and_write_for_one_worker(input_mat_file, taskid, start_row, end_row, total_column,
			mem_size_of_task, output_dir, job_id);

		start_row = end_row + 1;
	}

	t_end = time(NULL);
	t_elapse = difftime( t_end, t_start );
	printf("\t>> split_data total time - %d:%d:%d\n", (int)floor(t_elapse/3600.0), (int)floor(fmod(t_elapse,3600.0)/60.0), (int)fmod(t_elapse,60.0));

	// free space
	printf("free space\n"); fflush( stdout );
	free_VEC_INT( &part_sizes );
	free_VEC_INT2( &part_info_of_workers );

	printf ("Done.\n");
	return(EXIT_SUCCESS);
}


