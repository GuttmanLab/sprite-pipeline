/*
 * split_data_parallel.c
 *
 *  Created on: August 14, 2014
 *      Author: Wenyuan Li
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "matutils.h"
#include <mpi.h>

#define MAPI_USER_ABORT                         999

#define MASTER 0               /* taskid of first task */
#define PARTITION_TAG 1        /* message type: partition info */
#define ITERATION_BEGIN_TAG 2  /* message type: begin an iteration */
#define D_SAVE_RETURN_TAG 3    /* message type: return d_save info to master */
#define END_TAG 4              /* message type: end signal */
#define ERROR_TAG 99              /* message type: end signal */

void print_usage(FILE* stream)
{
	fprintf(stream, "Usage:\n");
	fprintf(stream, "\tsplit_data_parallel <input matrix file> <#rows/columns> <output directory> <#task> <memory size (MB) for each task> <job_id>\n");
	fprintf(stream, "\nVersion: 1.0\n");
	fprintf(stream, "\nAuthor: Wenyuan Li\n");
	fprintf(stream, "Date: August 14, 2014\n");
}

int main(int argc, char ** argv) {
	int	numtasks,              /* number of tasks in partition */
		taskid,                /* a task identifier */
		numworkers,            /* number of worker tasks */
		dest;                  /* task id of message destination */
	MPI_Status status;
	int i, ret, cnt_recv, finished;
	float x;
	time_t t_start, t_end;
	double t_elapse; // seconds
	char input_mat_file[MAXCHAR], temp_data_file[MAXCHAR];
	int total_column = atoi(argv[2]);
	float mem_size_of_task = (float)atof(argv[5]);
	char output_dir[MAXCHAR], job_id[MAXCHAR];
	VEC_INT part_sizes;
	int cnt_row, start_row, end_row, sub_start_row, sub_end_row;

	if (argc!=7) {
		fprintf(stderr,"Error: split_data_parallel arguments (%d) are insufficient.\n\n", argc);
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

    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    if (numtasks < 2 ) {
    		fprintf(stderr,"Error: Need at least two MPI tasks. Quitting...\n");
    		MPI_Abort(MPI_COMM_WORLD, MAPI_USER_ABORT);
    		exit(EXIT_SUCCESS);
	}
    numworkers = numtasks-1;

    /**************************** master task ************************************/
    if (taskid == MASTER) {
		printf("split_data_parallel arguments:\n");
		printf("\tinput file: %s\n", input_mat_file);
		printf("\ttotal_column: %d\n", total_column);
		printf("\toutput directory: %s\n", output_dir);
		printf("\t#task: %d\n", numtasks);
		printf("\tmemory size for each cpu: %g\n", mem_size_of_task);
		fflush( stdout );

		t_start = time(NULL);

		partitions2( total_column, numworkers, &part_sizes ); // #rows for each task
		printf("split_data_parallel has started with %d tasks (both master and %d workers).\n", numtasks, numworkers);
		printf("master: partition sizes (%d partitions): ", part_sizes.n);
		put_VEC_INT(stdout, part_sizes, ", ");
		fflush( stdout );
		if (numworkers!=part_sizes.n) {
			fprintf(stderr, "Error (master): #parts %d != numworkers %d", part_sizes.n, numworkers);
			MPI_Abort(MPI_COMM_WORLD, MAPI_USER_ABORT);
			exit(EXIT_FAILURE);
		}

		// send partition info to each worker
		printf("master: send partition info to each worker\n"); fflush( stdout );
		for (dest=1, start_row=1; dest<=numworkers; dest++) {
			// the i-th part
			// start_row and end_row are 1-based index
			i = dest - 1;
			end_row = start_row + part_sizes.v[i] - 1;
			ret = quick_check_partitions_by_mem_size(part_sizes.v[i], total_column, mem_size_of_task, &x);
			if (ret==FALSE) {
				fprintf(stderr, "Error (master): max_memory_size (%gMB) for each task cannot afford to compute on the matrix (size=%d rows X %d columns, each row as a partition needs %lf MB memory)\n", mem_size_of_task, part_sizes.v[i], total_column, x);
				MPI_Abort(MPI_COMM_WORLD, MAPI_USER_ABORT);
				exit(EXIT_FAILURE);
			}
			MPI_Send(&start_row, 1, MPI_INT, dest, PARTITION_TAG, MPI_COMM_WORLD);
			MPI_Send(&end_row, 1, MPI_INT, dest, PARTITION_TAG, MPI_COMM_WORLD);
			printf("master: worker %d - (%d, %d)\n", dest, start_row, end_row); fflush( stdout );
			start_row = end_row + 1;
		}

		// receive finish message from all workers
		printf("master: will receive finish message from all workers\n"); fflush( stdout );
		cnt_recv = 0;
		while (cnt_recv<numworkers) {
			// wait for data sent from workers
			MPI_Recv(&finished, 1, MPI_INT, MPI_ANY_SOURCE, END_TAG, MPI_COMM_WORLD, &status);
			dest =status.MPI_SOURCE;
			printf("\tmaster: receive finish message from worker %d\n", dest); fflush( stdout );
			cnt_recv++;
		}

		t_end = time(NULL);
		t_elapse = difftime( t_end, t_start );
		printf("\t>> split_data_parallel master computation total time - %d:%d:%d\n", (int)floor(t_elapse/3600.0), (int)floor(fmod(t_elapse,3600.0)/60.0), (int)fmod(t_elapse,60.0));

		// free space
		free_VEC_INT( &part_sizes );

		printf ("Done (master).\n"); fflush( stdout );
    }

    /**************************** worker task ************************************/
    if (taskid > MASTER) {
		// start_row and end_row are 1-based
		// receive info of which partition assigned to this worker task
		MPI_Recv(&start_row, 1, MPI_INT, MASTER, PARTITION_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&end_row, 1, MPI_INT, MASTER, PARTITION_TAG, MPI_COMM_WORLD, &status);
		cnt_row = end_row - start_row + 1;
		calc_partitions_by_mem_size(cnt_row, total_column, mem_size_of_task, &part_sizes);
		printf("worker %d/%d: %d - %d (%d parts), part sizes:", taskid,numworkers,start_row,end_row,part_sizes.n);
		put_VEC_INT(stdout, part_sizes, ", ");

		for (i=0, sub_start_row=start_row; i<part_sizes.n; i++) {
			// the i-th part
			// start_row and end_row are 1-based index
			sub_end_row = sub_start_row + part_sizes.v[i] - 1;
			// extract the partition and write it to the temporary file
			sprintf(temp_data_file,"%s.temp.%s.task%d.subpart%d",output_dir,job_id,taskid,i);
			printf("\tworker %d: read sub-part %d (%d, %d) and write it to file '%s'\n", taskid,
					i, sub_start_row, sub_end_row, temp_data_file); fflush( stdout );
			ret = read_write_part_of_MATRIX_FLOAT(input_mat_file, sub_start_row, sub_end_row, total_column, temp_data_file);
			sub_start_row = sub_end_row + 1;
		} // end of for "iterate on parts"

		// send finished message to Master
		finished = 1;
		MPI_Send(&finished, 1, MPI_INT, MASTER, END_TAG, MPI_COMM_WORLD);

		// free space
		free_VEC_INT( &part_sizes );

		printf ("Done (worker %d).\n", taskid);
    } // end of if "worker"

    MPI_Finalize();
	return(EXIT_SUCCESS);
}

