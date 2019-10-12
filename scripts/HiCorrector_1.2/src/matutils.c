/*
 * matutils.c
 *
 *  Created on: Jun 19, 2014
 *      Author: Wenyuan Li
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "matutils.h"

/*--------------------
 Routines
--------------------*/
/*---------------------------------
		Misc
-----------------------------------*/
void errorfile(char *fn, char *operation)
{
	char msg[MAXCHAR];
	sprintf(msg,"Error : Cannot %s file %s !\n",operation,fn);
	fprintf(stderr,msg);
	exit(-1);
}
void erroralloc(char *type, int count)
{
	char msg[MAXCHAR];
	sprintf(msg,"Error : Cannot allocate %d spaces for type %s !\n",count,type);
	fprintf(stderr,msg);
	exit(-1);
}
void remove_newlinechars( char* string )
{
	char D=0xD, A=0xA;
	char *p;
	/* DOS, Windows   0D 0A   (\r\n)
	   UNIX           0A      (\n)
	   Machintosh     0D      (\r)*/
	for ( p=string; *p!='\0'; p++ ) {
		if (*p==D) *p='\0';
		if (*p==A) *p='\0';
	}
}

/*---------------------------------
    Routines for VEC_INT
-----------------------------------*/
void init_VEC_INT(VEC_INT* vec)
{
	vec->n = 0;
	vec->v = NULL;
}
void free_VEC_INT(VEC_INT* vec)
{
	if (vec->v!=NULL) free( vec->v );
}
void create_VEC_INT(int n, VEC_INT* vec)
{
	vec->n = n;
	vec->v = (int*)malloc( n*sizeof(int) );
	if ( vec->v == NULL ) erroralloc("int",n);
}
void put_VEC_INT(FILE* stream, VEC_INT x, char* delim)
{
	int i;
	for (i=0; i<x.n-1; i++) fprintf( stream, "%d%s", x.v[i], delim );
	fprintf( stream, "%d\n", x.v[i] );
}
int max_VEC_INT(VEC_INT vec)
{
	int i, max=0;
	if (vec.v==NULL) return max;
	max = vec.v[0];
	for (i=0; i<vec.n; i++)
		if ( vec.v[i] > max ) max = vec.v[i];
	return max;
}


/*---------------------------------
    Routines for VEC_INT2
-----------------------------------*/
void init_VEC_INT2(VEC_INT2* vec)
{
	vec->n = 0;
	vec->v1 = NULL;
	vec->v2 = NULL;
}
void free_VEC_INT2(VEC_INT2* vec)
{
	if (vec->v1!=NULL) free( vec->v1 );
	if (vec->v2!=NULL) free( vec->v2 );
}
void create_VEC_INT2(int n, VEC_INT2* vec)
{
	vec->n = n;
	vec->v1 = (int*)malloc( n*sizeof(int) );
	if ( vec->v1 == NULL ) erroralloc("int",n);
	vec->v2 = (int*)malloc( n*sizeof(int) );
	if ( vec->v2 == NULL ) erroralloc("int",n);
}
void put_VEC_INT2(FILE* stream, VEC_INT2 x, char* delim)
{
	int i;
	for (i=0; i<x.n-1; i++) fprintf( stream, "(%d, %d)%s", x.v1[i], x.v2[i], delim );
	fprintf( stream, "(%d, %d)\n", x.v1[i], x.v2[i] );
}

/*---------------------------------
    Routines for VEC_FLOAT
-----------------------------------*/
void init_VEC_FLOAT(VEC_FLOAT* vec)
{
	vec->n = 0;
	vec->v = NULL;
}
void addnumber_VEC_FLOAT( float key, VEC_FLOAT* vec )
{
	int n;
	float *tmp;
	n = vec->n + 1;
	tmp = (float*)realloc( vec->v, n*sizeof(float) );
	if ( tmp == NULL ) erroralloc("float", n);
	vec->v = tmp;
	vec->v[n-1] = key;
	vec->n = n;
}
void free_VEC_FLOAT(VEC_FLOAT* vec)
{
	if (vec->v!=NULL) free( vec->v );
}
void create_VEC_FLOAT(int n, VEC_FLOAT* vec)
{
	vec->n = n;
	vec->v = (float*)malloc( n*sizeof(float) );
	if ( vec->v == NULL ) erroralloc("float",n);
}
void create_ones_VEC_FLOAT(int n, VEC_FLOAT* vec)
{
	int i;
	vec->n = n;
	vec->v = (float*)malloc( n*sizeof(float) );
	if ( vec->v == NULL ) erroralloc("float",n);
	for ( i=0; i<n; i++ ) vec->v[i] = 1;
}
int read_VEC_FLOAT(VEC_FLOAT* vec, char* file)
{
	FILE* stream;
	float v;
	int flag, i;
	init_VEC_FLOAT( vec );
	if( (stream = fopen( file, "r" )) == NULL ) { errorfile( file, "read" ); exit(EXIT_FAILURE); }
	i = 0;
	while( !feof( stream ) ) {
		flag = fscanf( stream, "%f", &v );
		if (flag>0) {
			addnumber_VEC_FLOAT(v, vec);
			i++;
		}
	}
	fclose(stream);
	return i;
}
int read_VEC_FLOAT_2(float* vec, int len, char* file)
{
	FILE* stream;
	float v;
	int flag, i;
	if( (stream = fopen( file, "r" )) == NULL ) { errorfile( file, "read" ); exit(EXIT_FAILURE); }
	i = 0;
	while( !feof( stream ) && i<len ) {
		flag = fscanf( stream, "%f", &v );
		if (flag>0) {
			vec[i] = v;
			i++;
		} else {
			fprintf(stderr, "Error 'read_VEC_FLOAT_2' read float wrong (i=%d)!\n",i);
		}
	}
	fclose(stream);
	return i;
}
void elementwise_multi_vector_VEC_FLOAT(VEC_FLOAT* dst, int start_idx,int end_idx,
		VEC_FLOAT src)
{
	int i, cnt_elements=end_idx-start_idx+1;
	if (cnt_elements!=src.n) {
		fprintf(stderr, "Error 'elementwise_multi_vector_VEC_FLOAT':\n\tlength of src (%d) doesn't match with (start_idx-end_idx+1)=%d\n", src.n, cnt_elements);
		exit(EXIT_FAILURE);
	}
	for (i=0; i<cnt_elements; i++)
		dst->v[i+start_idx-1] *= src.v[i];
}
// multiplication of elements of src1 from start_idx to end_idx with all elements of src2, then save
// to dst from start_idx to end_idx. It is supposed that src1 and dst have the same length.
// start_idx and end_idx are 1-based
void elementwise_multi_vector_VEC_FLOAT_2(float* src1, int start_idx, int end_idx,
		VEC_FLOAT src2, float* dst)
{
	int i, cnt_elements=end_idx-start_idx+1;
	if (cnt_elements!=src2.n) {
		fprintf(stderr, "Error 'elementwise_multi_vector_VEC_FLOAT':\n\tlength of src (%d) doesn't match with (start_idx-end_idx+1)=%d\n", src2.n, cnt_elements);
		exit(EXIT_FAILURE);
	}
	for (i=0; i<cnt_elements; i++)
		dst[i+start_idx-1] = src1[i+start_idx-1] * src2.v[i];
}

void assign_init_to_VEC_FLOAT(VEC_FLOAT* vec, float init_v)
{
	int i;
	for (i=0; i<vec->n; i++) vec->v[i] = init_v;
}

void write_VEC_FLOAT(VEC_FLOAT x, char* file)
{
	FILE* stream;
	int i;
	if( (stream = fopen( file, "w" )) == NULL ) { errorfile( file, "write" ); exit(EXIT_FAILURE); }
	for (i=0; i<x.n; i++) fprintf( stream, "%g\n", x.v[i] );
	//fprintf( stream, "\n");
	fclose(stream);
}

void write_VEC_FLOAT_2(float *x, int len, char* delim, char* file)
{
	FILE* stream;
	int i;
	if( (stream = fopen( file, "w" )) == NULL ) { errorfile( file, "write" ); exit(EXIT_FAILURE); }
	for (i=0; i<len-1; i++) fprintf( stream, "%g%s", x[i], delim );
	fprintf( stream, "%g\n", x[i] );
	fclose(stream);
}

void put_VEC_FLOAT(FILE* stream, VEC_FLOAT x, char* delim)
{
	int i;
	for (i=0; i<x.n-1; i++) fprintf( stream, "%g%s", x.v[i], delim );
	fprintf( stream, "%g\n", x.v[i] );
}

void put_VEC_FLOAT_2(FILE* stream, float* x, int len, char* delim)
{
	int i;
	for (i=0; i<len-1; i++) fprintf( stream, "%g%s", x[i], delim );
	fprintf( stream, "%g\n", x[i] );
}

// initialize a vector where all elements are ones.
void ones_VEC_FLOAT( VEC_FLOAT* vec )
{
	int i;
	for (i=0; i<vec->n; i++) vec->v[i]=1;
}

/* suppose dst is allocated space, before calling, index is 0-based
 * start_idx and start_idx are all indexes of dst vector*/
int copy_VEC_FLOAT(VEC_FLOAT *dst, int start_idx, int end_idx, VEC_FLOAT src)
{
	int i, j;
	if (src.n<(end_idx-start_idx+1)) {
		fprintf(stderr, "ERROR (copy_VEC_FLOAT): length (src.n=%d) is less than start=%d,end=%d!\n", src.n, start_idx, end_idx );
		return FALSE;
	}
	if (src.v==NULL) {
		fprintf(stderr, "ERROR (copy_VEC_FLOAT): src vector is NULL!\n");
		return FALSE;
	}
	for ( j=0, i=start_idx; i<=end_idx; i++, j++ ) dst->v[i] = src.v[j];
	return TRUE;
}

/*---------------------------------
    Routines for MATRIX_FLOAT
-----------------------------------*/
void init_MATRIX_FLOAT(MATRIX_FLOAT* mat)
{
	mat->start_row = 0;
	mat->end_row = 0;
	mat->cnt_row = 0;
	mat->real_cnt_row = 0;
	mat->total_column = 0;
	mat->v = NULL;
}

void create_MATRIX_FLOAT(int cnt_row, int total_column, MATRIX_FLOAT* mat)
{
	int i;
	mat->cnt_row = cnt_row;
	mat->real_cnt_row = cnt_row;
	mat->total_column = total_column;
	mat->v = (float**)malloc( cnt_row*sizeof(float*) );
	if ( mat->v == NULL ) erroralloc("float*",cnt_row);
	for ( i=0; i<cnt_row; i++ ) {
		mat->v[i] = (float*)malloc( total_column*sizeof(float) );
		if ( mat->v[i] == NULL ) erroralloc("float",total_column);
	}
}

void add_morerows_MATRIX_FLOAT(int n_additional_row, int total_column, MATRIX_FLOAT* mat)
{
	float **tmp;
	int i;
	tmp = (float**)realloc( mat->v, (mat->real_cnt_row+n_additional_row)*sizeof(float*) );
	if ( tmp == NULL ) erroralloc("float*",mat->real_cnt_row+n_additional_row);
	mat->v = tmp;
	for ( i=mat->real_cnt_row; i<mat->real_cnt_row+n_additional_row+n_additional_row; i++ ) {
		mat->v[i] = (float*)malloc( total_column*sizeof(float) );
		if ( mat->v[i] == NULL ) erroralloc("float",total_column);
	}
	mat->real_cnt_row += n_additional_row;
}

void create_ones_MATRIX_FLOAT(int cnt_row, int total_column, MATRIX_FLOAT* mat)
{
	int i,j;
	mat->cnt_row = cnt_row;
	mat->real_cnt_row= cnt_row;
	mat->total_column = total_column;
	mat->start_row = -1;
	mat->end_row = -1;
	mat->v = (float**)malloc( cnt_row*sizeof(float*) );
	if ( mat->v == NULL ) erroralloc("float*",cnt_row);
	for ( i=0; i<cnt_row; i++ ) {
		mat->v[i] = (float*)malloc( total_column*sizeof(float) );
		if ( mat->v[i] == NULL ) erroralloc("float",total_column);
		for ( j=0; j<total_column; j++ ) mat->v[i][j] = 1;
	}
}

void free_MATRIX_FLOAT( MATRIX_FLOAT* mat )
{
	int i;
	if (mat->v!=NULL) {
		for ( i=0; i<mat->real_cnt_row; i++ )
			if (mat->v[i]!=NULL) free( mat->v[i] );
		free( mat->v );
	}
}

//
// start_row and end_row are 1-based index
// read rows from "start_row" to "end_row"
// no header line in the input matrix file.
//
// start_row and int end_row are 1-based.
//
int read_MATRIX_FLOAT(char* file, int start_row, int end_row, int total_column,
		MATRIX_FLOAT* mat)
{
	FILE* stream;
	//char seps[] = " "; // a space
	char seps[] = "\t"; // a tab
	int cnt_row = end_row - start_row + 1;
	int i, j;
	char line[MAX_LINE_LENGTH], *token;
	free_MATRIX_FLOAT( mat );
	init_MATRIX_FLOAT( mat );
	create_MATRIX_FLOAT(cnt_row, total_column, mat);
	mat->start_row = start_row;
	mat->end_row = end_row;
	if( (stream = fopen( file, "r" )) == NULL ) {
		errorfile( file, "read" ); exit(0);
	}
	// skip the first (start_row-1) lines
	for (i=0; i<start_row-1; i++)
		fgets( line, MAX_LINE_LENGTH, stream );
	// read data starting from (start_row) line
	while( !feof( stream ) ) {
		for (i=0; i<cnt_row; i++) {
			if( fgets( line, MAX_LINE_LENGTH, stream ) != NULL) {
				remove_newlinechars( line );
				token = strtok( line, seps );
				j = 0;
				while( token != NULL && j<total_column ) {
					mat->v[i][j] = (float) atof(token);
					token = strtok( NULL, seps );
					j++;
				}
			} else {
				break;
			}
			if (j<total_column) {
				fprintf(stderr, "Warning: Line %d has %d elements which are not enough (%d lines)\n",
						start_row+i, j, total_column);
				return(EXIT_FAILURE);
			}
		}
		if (i<cnt_row) {
			fprintf(stderr, "Warning: Not enough lines %d which are not enough (%d lines)\n", i, cnt_row);
			return(EXIT_FAILURE);
		} else {
			break;
		}
	}
	fclose( stream );
	return(EXIT_SUCCESS);
}

//
// start_row and end_row are 1-based index
// read rows from "start_row" to "end_row"
// no header line in the input matrix file.
//
// start_row and int end_row are 1-based.
//
// if has_header_row==1, this file has 1st line which is header.
// if has_header_column==1, this file has 1st column is header.
int read_efficient_MATRIX_FLOAT_with_header(char* file, int start_row, int end_row, int total_column,
		int has_header_row, int has_header_column, MATRIX_FLOAT* mat)
{
	FILE* stream;
	//char seps[] = " "; // a space
	char seps[] = "\t"; // a tab
	int cnt_row = end_row - start_row + 1;
	int i, j;
	char line[MAX_LINE_LENGTH], *token;
	if (mat->v==NULL) {
		init_MATRIX_FLOAT( mat );
		create_MATRIX_FLOAT(cnt_row, total_column, mat);
	} else {
		if (mat->real_cnt_row<cnt_row) {
			add_morerows_MATRIX_FLOAT(cnt_row-mat->real_cnt_row, total_column, mat);
		}
		mat->cnt_row = cnt_row;
	}
	mat->start_row = start_row;
	mat->end_row = end_row;
	if( (stream = fopen( file, "r" )) == NULL ) {
		errorfile( file, "read" ); exit(0);
	}
	if (has_header_row==1)
		fgets( line, MAX_LINE_LENGTH, stream ); // skip header line
	// skip the first (start_row-1) lines
	for (i=0; i<start_row-1; i++)
		fgets( line, MAX_LINE_LENGTH, stream );
	// read data starting from (start_row) line
	while( !feof( stream ) ) {
		for (i=0; i<cnt_row; i++) {
			if( fgets( line, MAX_LINE_LENGTH, stream ) != NULL) {
				// read and process one line
				remove_newlinechars( line );
				token = strtok( line, seps );
				if (has_header_column==1)
					token = strtok( NULL, seps ); // skip 1st column which is header.
				j = 0;
				while( token != NULL && j<total_column ) {
					mat->v[i][j] = (float) atof(token);
					token = strtok( NULL, seps );
					j++;
				}
			} else {
				break;
			}
			if (j<total_column) {
				fprintf(stderr, "Warning: Line %d has %d elements which are not enough (%d lines)\n",
						start_row+i, j, total_column);
				return(EXIT_FAILURE);
			}
		}
		if (i<cnt_row) {
			fprintf(stderr, "Warning: Not enough lines %d which are not enough (%d lines)\n", i, cnt_row);
			return(EXIT_FAILURE);
		} else {
			break;
		}
	}
	fclose( stream );
	return(EXIT_SUCCESS);
}
//
// start_row and end_row are 1-based index
// read rows from "start_row" to "end_row"
// no header line in the input matrix file.
//
// start_row and int end_row are 1-based.
//
int read_efficient_MATRIX_FLOAT(char* file, int start_row, int end_row, int total_column,
		MATRIX_FLOAT* mat)
{
	FILE* stream;
	//char seps[] = " "; // a space
	char seps[] = "\t"; // a tab
	int cnt_row = end_row - start_row + 1;
	int i, j;
	char line[MAX_LINE_LENGTH], *token;
	if (mat->v==NULL) {
		init_MATRIX_FLOAT( mat );
		create_MATRIX_FLOAT(cnt_row, total_column, mat);
	} else {
		if (mat->real_cnt_row<cnt_row) {
			add_morerows_MATRIX_FLOAT(cnt_row-mat->real_cnt_row, total_column, mat);
		}
		mat->cnt_row = cnt_row;
	}
	mat->start_row = start_row;
	mat->end_row = end_row;
	if( (stream = fopen( file, "r" )) == NULL ) {
		errorfile( file, "read" ); exit(0);
	}
	// skip the first (start_row-1) lines
	for (i=0; i<start_row-1; i++)
		fgets( line, MAX_LINE_LENGTH, stream );
	// read data starting from (start_row) line
	while( !feof( stream ) ) {
		for (i=0; i<cnt_row; i++) {
			if( fgets( line, MAX_LINE_LENGTH, stream ) != NULL) {
				remove_newlinechars( line );
				token = strtok( line, seps );
				j = 0;
				while( token != NULL && j<total_column ) {
					mat->v[i][j] = (float) atof(token);
					token = strtok( NULL, seps );
					j++;
				}
			} else {
				break;
			}
			if (j<total_column) {
				fprintf(stderr, "Warning: Line %d has %d elements which are not enough (%d lines)\n",
						start_row+i, j, total_column);
				return(EXIT_FAILURE);
			}
		}
		if (i<cnt_row) {
			fprintf(stderr, "Warning: Not enough lines %d which are not enough (%d lines)\n", i, cnt_row);
			return(EXIT_FAILURE);
		} else {
			break;
		}
	}
	fclose( stream );
	return(EXIT_SUCCESS);
}

//
// start_row and end_row are 1-based index
// read rows from "start_row" to "end_row"
// no header line in the input matrix file.
//
// start_row and int end_row are 1-based.
//
int read_write_efficient_MATRIX_FLOAT(char* file, int start_row, int end_row, int total_column,
		MATRIX_FLOAT* mat, char* out_file)
{
	FILE *stream, *stream_out;
	//char seps[] = " "; // a space
	char seps[] = "\t"; // a tab
	int cnt_row = end_row - start_row + 1;
	int i, j;
	char line[MAX_LINE_LENGTH], *token;
	if (mat->v==NULL) {
		init_MATRIX_FLOAT( mat );
		create_MATRIX_FLOAT(cnt_row, total_column, mat);
	} else {
		if (mat->real_cnt_row<cnt_row) {
			add_morerows_MATRIX_FLOAT(cnt_row-mat->real_cnt_row, total_column, mat);
		}
		mat->cnt_row = cnt_row;
	}
	mat->start_row = start_row;
	mat->end_row = end_row;
	if( (stream = fopen( file, "r" )) == NULL ) {
		errorfile( file, "read" ); exit(0);
	}
	// skip the first (start_row-1) lines
	for (i=0; i<start_row-1; i++)
		fgets( line, MAX_LINE_LENGTH, stream );
	// read and write data starting from (start_row) line
	if( (stream_out = fopen( out_file, "w" )) == NULL ) {
		errorfile( file, "write" ); exit(0);
	}
	while( !feof( stream ) ) {
		for (i=0; i<cnt_row; i++) {
			if( fgets( line, MAX_LINE_LENGTH, stream ) != NULL) {
				fputs( line, stream_out );
				remove_newlinechars( line );
				token = strtok( line, seps );
				j = 0;
				while( token != NULL && j<total_column ) {
					mat->v[i][j] = (float) atof(token);
					token = strtok( NULL, seps );
					j++;
				}
			} else {
				break;
			}
			if (j<total_column) {
				fprintf(stderr, "Warning: Line %d has %d elements which are not enough (%d lines)\n",
						start_row+i, j, total_column);
				return(EXIT_FAILURE);
			}
		}
		if (i<cnt_row) {
			fprintf(stderr, "Warning: Not enough lines %d which are not enough (%d lines)\n", i, cnt_row);
			return(EXIT_FAILURE);
		} else {
			break;
		}
	}
	fclose( stream );
	fclose( stream_out );
	return(EXIT_SUCCESS);
}

//
// start_row and end_row are 1-based index
// read rows from "start_row" to "end_row"
// no header line in the input matrix file.
//
// start_row and int end_row are 1-based.
//
int read_write_part_of_MATRIX_FLOAT(char* file, int start_row, int end_row, int total_column, char* out_file)
{
	FILE *stream, *stream_out;
	int cnt_row = end_row - start_row + 1;
	int i;
	char line[MAX_LINE_LENGTH];
	if( (stream = fopen( file, "r" )) == NULL ) {
		errorfile( file, "read" ); exit(0);
	}
	// skip the first (start_row-1) lines
	for (i=0; i<start_row-1; i++)
		fgets( line, MAX_LINE_LENGTH, stream );
	// read and write data starting from (start_row) line
	if( (stream_out = fopen( out_file, "w" )) == NULL ) {
		errorfile( out_file, "write" ); exit(0);
	}
	while( !feof( stream ) ) {
		for (i=0; i<cnt_row; i++) {
			if( fgets( line, MAX_LINE_LENGTH, stream ) != NULL) {
				fputs( line, stream_out );
			} else {
				break;
			}
		}
		if (i<cnt_row) {
			fprintf(stderr, "Warning: Not enough lines %d which are not enough (%d lines)\n", i, cnt_row);
			return(EXIT_FAILURE);
		} else {
			break;
		}
	}
	fclose( stream );
	fclose( stream_out );
	return(EXIT_SUCCESS);
}

//
// no header line in the input matrix file.
// This function will read the cnt_row and then write to out_file
// It assumes "stream" is already open and it doesn't close it.
//
int read_write_part_of_MATRIX_FLOAT_alt(FILE* stream, int cnt_row, char* out_file)
{
	FILE *stream_out;
	int i;
	char line[MAX_LINE_LENGTH];
	// read and write data starting from (start_row) line
	if( (stream_out = fopen( out_file, "w" )) == NULL ) {
		errorfile( out_file, "write" ); exit(0);
	}
	while( !feof( stream ) ) {
		for (i=0; i<cnt_row; i++) {
			if( fgets( line, MAX_LINE_LENGTH, stream ) != NULL) {
				fputs( line, stream_out );
			} else {
				break;
			}
		}
		if (i<cnt_row) {
			fprintf(stderr, "Warning: Not enough lines %d which are not enough (%d lines)\n", i, cnt_row);
			return(EXIT_FAILURE);
		} else {
			break;
		}
	}
	fclose( stream_out );
	return(EXIT_SUCCESS);
}

void write_MATRIX_FLOAT(char* file, MATRIX_FLOAT mat)
{
	FILE* stream;
	int i,j;
	if( (stream = fopen( file, "w" )) == NULL ) { errorfile( file, "write" ); exit(0); }
	for (i=0; i<mat.cnt_row; i++) {
		for (j=0; j<mat.total_column-1; j++)
			fprintf(stream,"%f\t", mat.v[i][j]);
		fprintf(stream,"%f\n", mat.v[i][j]);
	}
	fclose(stream);
}

void append_MATRIX_FLOAT_text(char* file, MATRIX_FLOAT mat)
{
	FILE* stream;
	int i,j;
	if( (stream = fopen( file, "a" )) == NULL ) { errorfile( file, "append" ); exit(0); }
	for (i=0; i<mat.cnt_row; i++) {
		for (j=0; j<mat.total_column-1; j++)
			fprintf(stream,"%f\t", mat.v[i][j]);
		fprintf(stream,"%f\n", mat.v[i][j]);
	}
	fclose(stream);
}

void put_MATRIX_FLOAT(FILE* stream, MATRIX_FLOAT mat)
{
	int i,j;
	for (i=0; i<mat.cnt_row; i++) {
		for (j=0; j<mat.total_column-1; j++)
			fprintf(stream,"%.4f\t", mat.v[i][j]);
		fprintf(stream,"%.4f\n", mat.v[i][j]);
	}
}

float normp_float(float* v, int len, float p)
{
	int i;
	float ret=0;
	for (i=0; i<len; i++) ret+=(float)pow(v[i], p);
	return (float)pow(ret,1/p);
}

float norm1_float(float* v, int len)
{
	int i;
	float ret=0;
	for (i=0; i<len; i++) {
		ret+=v[i];
	}
	return ret;
}

float max_float(float* v, int len)
{
	int i;
	float max=v[0];
	for (i=0; i<len; i++)
		if (v[i]>max) max=v[i];
	return max;
}

void norm_rows_MATRIX_FLOAT(MATRIX_FLOAT mat, float p, VEC_FLOAT* norms)
{
	int i, cnt_zeros=0;
	if (norms->n!=mat.cnt_row) {
		free_VEC_FLOAT(norms);
		init_VEC_FLOAT(norms);
		create_VEC_FLOAT(mat.cnt_row, norms);
	}
	for (i=0; i<mat.cnt_row; i++) {
		if (p==-1)
			norms->v[i] = max_float(mat.v[i], mat.total_column);
		else if (p==1)
			norms->v[i] = norm1_float(mat.v[i], mat.total_column);
		else if (p>0)
			norms->v[i] = normp_float(mat.v[i], mat.total_column, p);
		else {
			fprintf(stderr, "Error: p<0 and p!=-1 is wrong!\n");
			exit(EXIT_FAILURE);
		}
		if (norms->v[i]==0) cnt_zeros++;
	}
	if (cnt_zeros!=0)
		printf("%d rows are all zeros.\n", cnt_zeros);
	fflush( stdout );
}

void elementwise_div_vector_MATRIX_FLOAT(MATRIX_FLOAT* mat, VEC_FLOAT norms)
{
	int i,j;
	if(mat->total_column!=norms.n) {
		fprintf(stderr,"Error 'elementwise_div_vector_MATRIX_FLOAT'\n\tLengths of matrix (%d) and vector (%d) doesn't match.\n", mat->total_column, norms.n);
		exit(EXIT_FAILURE);
	}
	for (i=0; i<mat->cnt_row; i++) {
		if (norms.v[i+mat->start_row-1]!=0) {
			for (j=0; j<mat->total_column; j++)
				if (norms.v[j]!=0) {
					mat->v[i][j] /= norms.v[i+mat->start_row-1]*norms.v[j];
				}
		}
	}
}

void elementwise_div_value_MATRIX_FLOAT(MATRIX_FLOAT* mat, float v)
{
	int i,j;
	if (v==0) return;
	for (i=0; i<mat->cnt_row; i++) {
		for (j=0; j<mat->total_column; j++)
			mat->v[i][j] /= v;
	}
}

void elementwise_div_vector_MATRIX_FLOAT_2(MATRIX_FLOAT* mat, float* norms, int norms_len)
{
	int i,j;
	if(mat->total_column!=norms_len) {
		fprintf(stderr,"Error 'elementwise_div_vector_MATRIX_FLOAT_2'\n\tLengths of matrix (%d) and vector (%d) doesn't match.\n", mat->total_column, norms_len);
		exit(EXIT_FAILURE);
	}
	for (i=0; i<mat->cnt_row; i++) {
		if (norms[i+mat->start_row-1]!=0) {
			for (j=0; j<mat->total_column; j++)
				if (norms[j]!=0) {
					mat->v[i][j] /= norms[i+mat->start_row-1]*norms[j];
				}
		}
	}
}

/*---------------------------------
    Routines for MATRIX_DOUBLE
-----------------------------------*/
void init_MATRIX_DOUBLE(MATRIX_DOUBLE* mat)
{
	mat->start_row = 0;
	mat->end_row = 0;
	mat->cnt_row = 0;
	mat->total_column = 0;
	mat->v = NULL;
}

void create_MATRIX_DOUBLE(int cnt_row, int total_column, MATRIX_DOUBLE* mat)
{
	int i,j;
	mat->cnt_row = cnt_row;
	mat->total_column = total_column;
	mat->v = (double**)malloc( cnt_row*sizeof(double*) );
	if ( mat->v == NULL ) erroralloc("double*",cnt_row);
	for ( i=0; i<cnt_row; i++ ) {
		mat->v[i] = (double*)malloc( total_column*sizeof(double) );
		if ( mat->v[i] == NULL ) erroralloc("double",total_column);
		for ( j=0; j<total_column; j++ ) mat->v[i][j] = 0;
	}
}

void free_MATRIX_DOUBLE( MATRIX_DOUBLE* mat )
{
	int i;
	if (mat->v!=NULL) {
		for ( i=0; i<mat->cnt_row; i++ )
			if (mat->v[i]!=NULL) free( mat->v[i] );
		free( mat->v );
	}
}

/*---------------------------------
    Routines for calculation of partitions
-----------------------------------*/
// Given a fixed partition size, get size of each partition.
// So the number of partitions may vary.
void partitions( int n, int part_size, VEC_INT* part_sizes )
{
	int n_part =  (int)floor(n/(float)part_size);
	int last_part_size = n - n_part * part_size;
	int i, last_part_additional_size, avg_more_size, remainder;
	if (last_part_size>0) {
		create_VEC_INT( n_part+1, part_sizes);
		for (i=0; i<n_part; i++)
			part_sizes->v[i] = part_size;
		part_sizes->v[i] = last_part_size;
		n_part++;
	} else {
		create_VEC_INT( n_part, part_sizes);
		for (i=0; i<n_part; i++)
			part_sizes->v[i] = part_size;
	}
	last_part_additional_size = part_sizes->v[n_part-1]-part_size;
	if (last_part_additional_size>1) {
		// the following can make part_size distribution is smoother.
		if (n_part>1) {
			avg_more_size =(int)floor(last_part_additional_size/(float)(n_part-1));
			remainder = last_part_additional_size%(n_part-1);
			if (avg_more_size==0) {
				for (i=0; i<remainder; i++)
					part_sizes->v[i]++;
				part_sizes->v[n_part-1] -= remainder;
			} else {
				for (i=0; i<n_part-1; i++)
					part_sizes->v[i]++;
				part_sizes->v[n_part-1] -= n_part-1;

				// do this process again, if last_part_additional_size is still large.
				// it is impossible that avg_more_size>0 anymore this time.
				last_part_additional_size = part_sizes->v[n_part-1]-part_size;
				if (last_part_additional_size>1) {
					avg_more_size =(int)floor(last_part_additional_size/(float)(n_part-1));
					remainder = last_part_additional_size/(n_part-1);
					if (avg_more_size==0) {
						for (i=0; i<remainder; i++)
							part_sizes->v[i]++;
						part_sizes->v[n_part-1] -= remainder;
					}
				}
			}
		}
	}
}
// Given a fixed number of partitions (n_part), get size of each
// partition. So the size of last partition may vary.
void partitions2( int n, int n_part, VEC_INT* part_sizes )
{
	int part_size = (int)floor(n/(float)n_part);
	int last_part_additional_size = n - n_part * part_size;
	int i, avg_more_size, remainder;
	create_VEC_INT( n_part, part_sizes);
	for (i=0; i<n_part-1; i++)
		part_sizes->v[i] = part_size;
	part_sizes->v[i] = part_size + last_part_additional_size;
	if (last_part_additional_size>1) {
		// the following can make part_size distribution is smoother.
		if (n_part>1) {
			avg_more_size =(int)floor(last_part_additional_size/(float)(n_part-1));
			remainder = last_part_additional_size%(n_part-1);
			if (avg_more_size==0) {
				for (i=0; i<remainder; i++)
					part_sizes->v[i]++;
				part_sizes->v[n_part-1] -= remainder;
			} else {
				for (i=0; i<n_part-1; i++)
					part_sizes->v[i]++;
				part_sizes->v[n_part-1] -= n_part-1;

				// do this process again, if the last_part_additional_size is still large.
				// it is impossible that avg_more_size>0 anymore this time.
				last_part_additional_size = part_sizes->v[n_part-1]-part_size;
				if (last_part_additional_size>1) {
					avg_more_size =(int)floor(last_part_additional_size/(float)(n_part-1));
					remainder = last_part_additional_size/(n_part-1);
					if (avg_more_size==0) {
						for (i=0; i<remainder; i++)
							part_sizes->v[i]++;
						part_sizes->v[n_part-1] -= remainder;
					}
				}
			}
		}
	}
}
// float is always 4 bytes in both 32bit and 64bit machine: see http://docs.oracle.com/cd/E18752_01/html/817-6223/chp-typeopexpr-2.html
//
// x is the part_size (i.e., number of rows). Following are bytes:
//
// line string used (each number occupies 25 char) when reading data: 25*total_column
// matrix loaded: 4*x*total_column
// d_save (2 vectors, one for used, the other for backup): 4*2*total_column
// d_temp (norms of rows): 4*x
//
// Memory used by data is (bytes): total_column*(4*x+33) + 4*x (bytes)
// Memory used by codes supposed to be 10MB: 50*1024*1024 (bytes)
//
// Given a max memory size that are available (T MB), we get the part_size (number of rows):
// x=(1024*1024*(T-1)-33*total_column)/(4*total_column+1);
//
// Input:
//   total_mem_size: in MB
void calc_partitions_by_mem_size(int cnt_rows, int total_column, float max_mem_size,
		VEC_INT* part_sizes)
{
	int part_size;
	double x;
	x = (1024.0*1024.0*(max_mem_size-MEM_NOT_CONSIDERED)-33.0*total_column)/(4.0*total_column+1.0);
	if (x<1.0) {
		x = 1.0 + ((4.0*total_column+1.0)*1 +33.0*total_column)/1024.0/1024.0;
		fprintf(stderr, "Error: max_memory_size (%gMB) cannot afford to compute on the matrix (size=%d rows X %d columns, each row as a partition needs %lf MB memory)\n", max_mem_size, cnt_rows, total_column, x);
		exit(EXIT_FAILURE);
	}
	part_size = floor(x);
	partitions( cnt_rows, part_size, part_sizes );
}
int quick_check_partitions_by_mem_size(int cnt_rows, int total_column, float max_mem_size, float *x)
{
	*x = (1024.0*1024.0*(max_mem_size-MEM_NOT_CONSIDERED)-33.0*total_column)/(4.0*total_column+1.0);
	if (*x<1) {
		*x = 1.0 + ((4.0*total_column+1.0)*1 +33.0*total_column)/1024.0/1024.0;
		return FALSE;
	} else {
		*x = 1.0 + ((4.0*total_column+1.0)*1 +33.0*total_column)/1024.0/1024.0;
		return TRUE;
	}
}
