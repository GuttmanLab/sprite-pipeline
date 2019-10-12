/*
 * ic_mep.c
 *
 *  Created on: Jun 26, 2014
 *  Revised on: Oct 17, 2014
 *		1. Add the function that can read input matrix file with header line or column

 *      Author: Wenyuan Li
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "matutils.h"
#include "my_getopt-1.5/getopt.h"

#define MAPI_USER_ABORT                         999

#define MASTER 0               /* taskid of first task */
#define PARTITION_TAG 1        /* message type: partition info */
#define ITERATION_BEGIN_TAG 2  /* message type: begin an iteration */
#define D_SAVE_RETURN_TAG 3    /* message type: return d_save info to master */
#define END_TAG 4              /* message type: end signal */
#define ERROR_TAG 99              /* message type: end signal */

#define DEFAULT_HAS_HEADER_LINE 0 // default value for having header line in input matrix file
#define DEFAULT_HAS_HEADER_COLUMN 0 // default value for having header column in input matrix file
#define DEFAULT_MAX_ITERATION 10 // default value for maximum number of iterations

/*-----------------------------------
command line parser's data structure
-----------------------------------*/
struct option longopts[] =
{
	/* name,                      has_arg,          flag, val */ /* longind */
	{ "help",                     no_argument,       0,   'h' }, /*       0 */
	{ "inputFile",                required_argument, 0,   'I' }, /*       1 */
	{ "useSplitInputFiles",       no_argument,       0,   's' }, /*       2 */
	{ "numRows",                  required_argument, 0,   'N' }, /*       3 */
	{ "numTask",                  required_argument, 0,   'T' }, /*       4 */
	{ "memSizePerTask",           required_argument, 0,   'M' }, /*       5 */
	{ "maxIteration",             required_argument, 0,   'r' }, /*       6 */
	{ "jobID",                    required_argument, 0,   'J' }, /*       7 */
	{ "outputFile",               required_argument, 0,   'O' }, /*       8 */
	{ "hasHeaderRow",             required_argument, 0,   'z' }, /*       3 */
	{ "hasHeaderColumn",          required_argument, 0,   'Z' }, /*       3 */
	/* end-of-list marker */
	{ 0, 0, 0, 0 }
};
char *shortopts = "hsI:N:T:M:r:J:O:z:Z:"; /* short options string */
int longind = 0; /* long option list index */

void print_usage(FILE* stream)
{
	fprintf( stream, "\nUsage: ic_mep [options] ...\n"
		             "Options:\n"
					 "\t-h or --help          :show this message and exit\n"
					 "\t--inputFile=FILE      :input contact matrix file name. Not required when using split files.\n"
					 "\t--hasHeaderRow=NUM    :input matrix file has header line or not. 1 means 'has'; 0 otherwise. default (%d).\n"
					 "\t--hasHeaderColumn=NUM :input matrix file has header column or not. 1 means 'has'; 0 otherwise. default (%d).\n"
					 "\t--useSplitInputFiles  :Use split contact matrix files.\n"
					 "\t--outputFile=FILE     :output file for bias vector, REQUIRED.\n"
					 "\t--numRows=NUM         :number of rows in the contact matrix, REQUIRED.\n"
					 "\t--numTask=NUM         :number of all tasks (including manager task), REQUIRED.\n"
					 "\t--memSizePerTask=NUM  :size (mega bytes) of memory used for each task, REQUIRED.\n"
					 "\t--jobID=STRING        :job identity string for giving temporary file names, REQUIRED.\n"
					 "\t--maxIteration=NUM    :maximum number of iterations, default (%d).\n",
			DEFAULT_HAS_HEADER_LINE, DEFAULT_HAS_HEADER_COLUMN, DEFAULT_MAX_ITERATION);
		
	fprintf(stream, "\nVersion: 1.1\n");
	fprintf(stream, "Author: Wenyuan Li\n");
	fprintf(stream, "Date: October 17, 2014\n");
	fprintf( stream, "Please email XJZHOU@usc.edu or WEL@usc.edu\n");
	fprintf( stream, "Further refers to http:\/\/zhoulab.usc.edu for help.\n\n");
}

int main(int argc, char ** argv) {
	int opt; /* during argument parsing, opt contains the return value from getopt() */
	int longind = 0; /* long option list index */
	int	numtasks=-1,              /* number of tasks in partition */
		taskid,                /* a task identifier */
		numworkers,            /* number of worker tasks */
		dest;                  /* task id of message destination */
	MPI_Status status;
	float p = 1.0;
	int useSplitInputFiles = 0; // default is not to use split input files
	int total_column = -1;
	int has_header_line = DEFAULT_HAS_HEADER_LINE;
	int has_header_column = DEFAULT_HAS_HEADER_COLUMN;
	int maxiter = DEFAULT_MAX_ITERATION;  // default of iterations is 10
	float mem_size_of_task = -1;
	char job_id[MAXCHAR]="";
	char input_mat_file[MAXCHAR]="";
	char output_file[MAXCHAR]="";
	int i, ret, part_size, cnt_recv, size;
	float x;
	time_t t_start0, t_start1, t_end1, t_end2;
	double t_elapse;
	char temp_file[MAXCHAR], output_dir[MAXCHAR], temp_data_file[MAXCHAR];
	char *pt;
	VEC_INT part_sizes;
	int cnt_row, start_row, end_row, sub_start_row, sub_end_row, iter, iter_msg;
	MATRIX_FLOAT mat; // submatrix loaded, its size is "4*part_size*total_column" bytes (float is 4 bytes in both 32bit and 64bit machine: see http://docs.oracle.com/cd/E18752_01/html/817-6223/chp-typeopexpr-2.html)
	MATRIX_FLOAT d_save; // accumulated bias factors (2 vectors, one for used, the other for backup), its size is"4*2*total_column" bytes
	VEC_FLOAT d_temp, d_save_master; // norms of rows, its size is "4*part_size" bytes
	float *d_save_use_ptr, *d_save_new_ptr;
	VEC_INT2 part_info_of_workers;

	/* parse all options from the command line */
	while ((opt=getopt_long(argc, argv, shortopts, longopts, &longind)) != -1) {
		switch (opt) {
			case 'h': /* --help */
				print_usage( stdout );
				exit(EXIT_FAILURE);
			case 's': /* --useSplitInputFiles */
				useSplitInputFiles = 1;
				break;
			case 'I': /* --inputFile */
				strcpy( input_mat_file, optarg );
				break;
			case 'N': /* --numRows */
				total_column = atoi( optarg );
				break;
			case 'T': /* --numTask */
				numtasks = atoi( optarg );
				break;
			case 'M': /* --memSizePerTask */
				mem_size_of_task = (float)atof( optarg );
				break;
			case 'r': /* --maxIteration */
				maxiter = atoi( optarg );
				break;
			case 'J': /* --jobID */
				strcpy( job_id, optarg );
				break;
			case 'O': /* --outputFile */
				strcpy( output_file, optarg );
				break;
			case 'z': /* --hasHeaderRow */
				has_header_line = atoi( optarg );
				break;
			case 'Z': /* --hasHeaderColumn */
				has_header_column = atoi( optarg );
				break;
			default: /* something unexpected has happened */
				print_usage( stderr );
				exit( EXIT_FAILURE );
		}
	}
	if (useSplitInputFiles==0 && strlen(input_mat_file)==0) {
		fprintf(stderr, "\nError: --inputFile=FILE is required!\n");
		print_usage( stderr );
		exit( EXIT_FAILURE );
	}
	if (strlen(output_file)==0) {
		fprintf(stderr, "\nError: --outputFile=FILE is required!\n");
		print_usage( stderr );
		exit( EXIT_FAILURE );
	}
	if (total_column==-1) {
		fprintf(stderr, "\nError: --numRows=NUM is required!\n");
		print_usage( stderr );
		exit( EXIT_FAILURE );
	}
	if (numtasks==-1) {
		fprintf(stderr, "\nError: --numTask=NUM is required!\n");
		print_usage( stderr );
		exit( EXIT_FAILURE );
	}
    if (numtasks < 2 ) {
    		fprintf(stderr,"Error: Need at least two MPI tasks!\n");
    		exit(EXIT_SUCCESS);
	}
	if (mem_size_of_task==-1) {
		fprintf(stderr, "\nError: --memSizePerTask=NUM is required!\n");
		print_usage( stderr );
		exit( EXIT_FAILURE );
	}
	if (strlen(job_id)==0) {
		fprintf(stderr, "\nError: --jobID=STRING is required!\n");
		print_usage( stderr );
		exit( EXIT_FAILURE );
	}
	strcpy(output_dir, output_file);
	pt = strrchr(output_dir, '/');
	if (pt==NULL)
		strcpy(output_dir, "./");
	else {
		*pt = '\0';
		strcat(output_dir,"/");
	}

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    numworkers = numtasks-1;

    /**************************** master task ************************************/
    if (taskid == MASTER) {
		printf("ic_mep arguments:\n");
		if ( useSplitInputFiles == 1 )
			printf("\tuse split input contact matrix file\n");
		else {
			printf("\tinput file: %s\n", input_mat_file);
			printf("\t\thas header line: %d\n", has_header_line);
			printf("\t\thas header column: %d\n", has_header_column);
		}
		printf("\toutput file: %s\n", output_file);
		printf("\tnumber of rows: %d\n", total_column);
		printf("\t#task: %d\n", numtasks);
		printf("\tmemory size for each cpu: %g MB\n", mem_size_of_task);
		printf("\tmaximum number of iterations: %d\n", maxiter);
		fflush( stdout );

		t_start0 = time(NULL);

		partitions2( total_column, numworkers, &part_sizes ); // #rows for each task
		printf("ic_mep has started with %d tasks (both master and %d workers).\n", numtasks, numworkers);
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
		create_VEC_INT2(numworkers, &part_info_of_workers);
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
			part_info_of_workers.v1[dest-1] = start_row;
			part_info_of_workers.v2[dest-1] = end_row;
			printf("master: worker %d - (%d, %d)\n", dest, start_row, end_row); fflush( stdout );
			start_row = end_row + 1;
		}

		// begin to compute (iterations)
		printf("master: begin to compute\n"); fflush( stdout );
		part_size = max_VEC_INT(part_sizes);
		create_VEC_FLOAT( part_size, &d_temp );
		create_ones_VEC_FLOAT( total_column, &d_save_master );
		for (iter=0; iter<maxiter; iter++) {
			printf("+++++++++ iter %d +++++++++\n", iter); fflush( stdout );
			t_start1 = time(NULL);
			// notify each worker to begin an iteration of computation
			for (dest=1; dest<=numworkers; dest++) {
				printf("\tmaster: iter %d - step 1: notify worker %d to start\n", iter, dest); fflush( stdout );
				MPI_Send(&iter, 1, MPI_INT, dest, ITERATION_BEGIN_TAG, MPI_COMM_WORLD);
			}

			// receive d_save from all workers
			printf("\tmaster: iter %d - step 2: will receive data from all workers\n", iter); fflush( stdout );
			cnt_recv = 0;
			while (cnt_recv<numworkers) {
				// wait for data sent from workers
				MPI_Recv(&d_temp.v[0], part_size, MPI_FLOAT, MPI_ANY_SOURCE, D_SAVE_RETURN_TAG, MPI_COMM_WORLD, &status);
				dest =status.MPI_SOURCE;
				MPI_Get_count( &status, MPI_FLOAT, &size );
				if (size!=(part_info_of_workers.v2[dest-1]-part_info_of_workers.v1[dest-1]+1)) {
					fprintf(stderr, "Error (master): received data length doesn't agree with expected!\n");
					fflush( stderr );
					MPI_Abort(MPI_COMM_WORLD, MAPI_USER_ABORT);
					exit(EXIT_FAILURE);
				}
				ret = copy_VEC_FLOAT(&d_save_master, part_info_of_workers.v1[dest-1]-1,
						part_info_of_workers.v2[dest-1]-1, d_temp);
				if (ret==FALSE) {
					fprintf(stderr, "Error (master): copy_VEC_FLOAT, Length (from worker %d) doesn't agree or null vector!\n",
							dest);
					MPI_Abort(MPI_COMM_WORLD, MAPI_USER_ABORT);
					exit(EXIT_FAILURE);
				}
				printf("\t\tmaster: iter %d - receive %d float from worker %d\n", iter, size, dest);
				fflush( stdout );
				cnt_recv++;
			}
			// write "d_save_master" to file
			printf("\tmaster: iter %d - step 3: write collected data to temp file\n", iter); fflush( stdout );
			sprintf(temp_file,"%s.temp.%s",output_dir,job_id);
			// sprintf(temp_file,"%s.temp.%s.%d",output_dir,job_id,iter);
			write_VEC_FLOAT(d_save_master, temp_file);

			t_end1 = time(NULL);
			t_elapse = difftime( t_end1, t_start1 );
			printf("\tmaster: iter %d spent time: %d:%d:%d\n", iter, (int)floor(t_elapse/3600.0), (int)floor(fmod(t_elapse,3600.0)/60.0), (int)fmod(t_elapse,60.0));
		}

		t_end2 = time(NULL);
		t_elapse = difftime( t_end2, t_start0 );
		printf("\t>> ic_mep master computation total time - %d:%d:%d\n", (int)floor(t_elapse/3600.0), (int)floor(fmod(t_elapse,3600.0)/60.0), (int)fmod(t_elapse,60.0));

		// write results to output file
		t_start1 = time(NULL);
		printf("master: write results to output file '%s'\n", output_file); fflush( stdout );
		write_VEC_FLOAT(d_save_master, output_file);
		t_end2 = time(NULL);
		t_elapse = difftime( t_end2, t_start1 );
		printf("\tmaster: write results I/O time: %d:%d:%d\n", (int)floor(t_elapse/3600.0), (int)floor(fmod(t_elapse,3600.0)/60.0), (int)fmod(t_elapse,60.0));

		// free space
		printf("master: free space\n"); fflush( stdout );
		free_VEC_INT( &part_sizes );
		free_VEC_INT2( &part_info_of_workers );
		free_VEC_FLOAT( &d_temp );
		free_VEC_FLOAT( &d_save_master );

		printf ("Done (master).\n");
    }

    /**************************** worker task ************************************/
    if (taskid > MASTER) {
		// start_row and end_row are 1-based
		// receive info of which partition assigned to this worker task
		MPI_Recv(&start_row, 1, MPI_INT, MASTER, PARTITION_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&end_row, 1, MPI_INT, MASTER, PARTITION_TAG, MPI_COMM_WORLD, &status);
		cnt_row = end_row - start_row + 1;
		calc_partitions_by_mem_size(cnt_row, total_column, mem_size_of_task, &part_sizes);
		part_size = max_VEC_INT(part_sizes);
		printf("worker %d/%d: %d - %d (%d parts), part sizes:", taskid,numworkers,start_row,end_row,part_sizes.n);
		put_VEC_INT(stdout, part_sizes, ", ");

		init_MATRIX_FLOAT( &mat );
		create_ones_VEC_FLOAT( part_size, &d_temp );
		create_ones_MATRIX_FLOAT( 2, total_column, &d_save);

		d_save_use_ptr = d_save.v[0];
		d_save_new_ptr = d_save.v[1];

		for (iter=0; iter<maxiter; iter++) {
			// Wait to receive ITERATION_BEGIN message from master
			MPI_Recv(&iter_msg, 1, MPI_INT, MASTER, ITERATION_BEGIN_TAG, MPI_COMM_WORLD, &status);
			printf("\tworker %d: iter %d - ITERATION_BEGIN message is received from master\n", taskid, iter); fflush( stdout );
			if(iter!=iter_msg) {
				fprintf(stderr, "Error (worker %d): Iteration %d doesn't synchronize with Iteration %d from master!\n",
						taskid, iter, iter_msg);
				fflush( stderr );
				MPI_Abort(MPI_COMM_WORLD, MAPI_USER_ABORT);
				exit(EXIT_FAILURE);
			}

			// read d_save from temporary file
			if (iter!=0) {
				sprintf(temp_file,"%s.temp.%s",output_dir,job_id);
				printf("\tworker %d: iter %d - read bias vector from temp file (%s)\n", taskid, iter, temp_file); fflush( stdout );
				ret = read_VEC_FLOAT_2(d_save_use_ptr, total_column, temp_file);
			}

			for (i=0, sub_start_row=start_row; i<part_sizes.n; i++) {
				// the i-th part
				// start_row and end_row are 1-based index
				sub_end_row = sub_start_row + part_sizes.v[i] - 1;
				if ( useSplitInputFiles == 0 ) {
					// read single file of contact matrix
					printf("\tworker %d: iter %d - read '%s' and compute sub-part %d (%d, %d)\n", taskid,
							iter, input_mat_file, i, sub_start_row, sub_end_row); fflush( stdout );
					ret = read_efficient_MATRIX_FLOAT_with_header(input_mat_file, sub_start_row, sub_end_row, total_column, has_header_line, has_header_column, &mat);
				} else {
					// read split files of contact matrix
					sprintf(temp_data_file,"%s.temp.%s.task%d.subpart%d",output_dir,job_id,taskid,i);
					printf("\tworker %d: iter %d - read data sub-part file '%s', and compute sub-part %d (%d, %d)\n", taskid,
							iter, temp_data_file, i, sub_start_row, sub_end_row); fflush( stdout );
					ret = read_efficient_MATRIX_FLOAT(temp_data_file, 1, part_sizes.v[i], total_column, &mat);
					mat.start_row = sub_start_row;
					mat.end_row = sub_end_row;
					mat.cnt_row = part_sizes.v[i];
				}
				if (ret==EXIT_SUCCESS) {
					// update matrix: B = diag(1./d_save)*A*diag(1./d_save);
					if (iter!=0)
						elementwise_div_vector_MATRIX_FLOAT_2( &mat, d_save_use_ptr, total_column);
					// get new norms from updated matrix: d_temp = normp(B,p);
					norm_rows_MATRIX_FLOAT(mat, p, &d_temp);
					// update accumulated factors: d_save = d_save.*d_temp;
					elementwise_multi_vector_VEC_FLOAT_2(d_save_use_ptr, sub_start_row, sub_end_row, d_temp, d_save_new_ptr);
				} else {
					fprintf(stderr, "Error (worker %d): Wrong reading '%s' on lines %d - %d\n",
							taskid, input_mat_file, sub_start_row, sub_end_row);
					fflush( stderr );
					MPI_Abort(MPI_COMM_WORLD, MAPI_USER_ABORT);
					exit(EXIT_FAILURE);
				}
				sub_start_row = sub_end_row + 1;
			} // end of for "iterate on parts"

			// send updated d_save to Master for collection
			printf("\tworker %d: iter %d - send updated data to master\n", taskid, iter); fflush( stdout );
			MPI_Send(d_save_new_ptr+start_row-1, cnt_row, MPI_FLOAT, MASTER, D_SAVE_RETURN_TAG, MPI_COMM_WORLD);
		} // end of for "iterate on iterations"

		// free space
		printf("worker %d: free space\n", taskid); fflush( stdout );
		free_VEC_INT( &part_sizes );
		free_MATRIX_FLOAT( &mat );
		free_MATRIX_FLOAT( &d_save );
		free_VEC_FLOAT( &d_temp );

		printf ("Done (worker %d).\n", taskid);
    } // end of if "worker"

    MPI_Finalize();
	return(EXIT_SUCCESS);
}

