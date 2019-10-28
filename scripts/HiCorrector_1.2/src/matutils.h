/*
 * matutils.h
 *
 *  Created on: Jun 19, 2014
 *      Author: Wenyuan Li
 */

#ifndef MATUTILS_H_
#define MATUTILS_H_

/*---------------------------
  Data Types
---------------------------*/
typedef int BOOL;
#define TRUE 1
#define FALSE 0

/*---------------------------
  Constants
---------------------------*/
#define MAXCHAR 1024
#define MAX_LINE_LENGTH 1000000
#define MAXCHAR_ID 64
#define MEM_NOT_CONSIDERED 0  // in MB

/*---------------------------
  Data structure
---------------------------*/
typedef struct _vector_i {
	int n;
	int* v;
} VEC_INT;

typedef struct _vector_i2 {
	int n;
	int* v1;
	int* v2;
} VEC_INT2;

typedef struct _vector_f {
	int n;
	float* v;
} VEC_FLOAT;

typedef struct _vector_d {
	int n;
	double* v;
} VEC_DOUBLE;

typedef struct _matrix_f {
	int start_row;
	int end_row;
	int cnt_row;
	int total_column;
	int real_cnt_row;
	float** v;
} MATRIX_FLOAT;

typedef struct _matrix_d {
	int start_row;
	int end_row;
	int cnt_row;
	int total_column;
	double** v;
} MATRIX_DOUBLE;



/*---------------------------
  Routines
---------------------------*/
int getline_gnu(char **lineptr, size_t *n, FILE *stream);

void init_VEC_INT(VEC_INT* vec);
void free_VEC_INT(VEC_INT* vec);
void create_VEC_INT(int n, VEC_INT* vec);
void put_VEC_INT(FILE* stream, VEC_INT x, char* delim);
int max_VEC_INT(VEC_INT vec);

void init_VEC_INT2(VEC_INT2* vec);
void free_VEC_INT2(VEC_INT2* vec);
void create_VEC_INT2(int n, VEC_INT2* vec);
void put_VEC_INT2(FILE* stream, VEC_INT2 x, char* delim);

void init_VEC_FLOAT(VEC_FLOAT* vec);
void create_VEC_FLOAT(int n, VEC_FLOAT* vec);
void create_ones_VEC_FLOAT(int n, VEC_FLOAT* vec);
void addnumber_VEC_FLOAT( float key, VEC_FLOAT* vec );
void free_VEC_FLOAT(VEC_FLOAT* vec);
int read_VEC_FLOAT(VEC_FLOAT* vec, char* file);
int read_VEC_FLOAT_2(float* vec, int len, char* file);
void assign_init_to_VEC_FLOAT(VEC_FLOAT* vec, float init_v);
void write_VEC_FLOAT(VEC_FLOAT x, char* file);
void write_VEC_FLOAT_2(float *x, int len, char* delim, char* file);
void put_VEC_FLOAT(FILE* stream, VEC_FLOAT x, char* delim);
void put_VEC_FLOAT_2(FILE* stream, float* x, int len, char* delim);
void ones_VEC_FLOAT( VEC_FLOAT* vec );
void elementwise_div_vector_MATRIX_FLOAT(MATRIX_FLOAT* mat, VEC_FLOAT norms);
void elementwise_div_value_MATRIX_FLOAT(MATRIX_FLOAT* mat, float v);
void elementwise_div_vector_MATRIX_FLOAT_2(MATRIX_FLOAT* mat, float* norms, int norms_len);
void elementwise_multi_vector_VEC_FLOAT(VEC_FLOAT* dst, int start_idx,int end_idx,
		VEC_FLOAT src);
void elementwise_multi_vector_VEC_FLOAT_2(float* src1, int start_idx,int end_idx,
		VEC_FLOAT src2, float* dst);
/* suppose dst is allocated space, before calling, index is 0-based
 * start_idx and start_idx are all indexes of dst vector*/
int copy_VEC_FLOAT(VEC_FLOAT *dst, int start_idx, int end_idx, VEC_FLOAT src);

void init_MATRIX_FLOAT(MATRIX_FLOAT* mat);
void create_ones_MATRIX_FLOAT(int cnt_row, int total_column, MATRIX_FLOAT* mat);
void add_morerows_MATRIX_FLOAT(int n_additional_row, int total_column, MATRIX_FLOAT* mat);
void free_MATRIX_FLOAT( MATRIX_FLOAT* mat );
int read_MATRIX_FLOAT(char* file, int start_row, int end_row, int total_column,
		MATRIX_FLOAT* mat);
int read_efficient_MATRIX_FLOAT(char* file, int start_row, int end_row, int total_column,
		MATRIX_FLOAT* mat);
int read_efficient_MATRIX_FLOAT_with_header(char* file, int start_row, int end_row, int total_column,
		int has_header_row, int has_header_column, MATRIX_FLOAT* mat);
int read_write_efficient_MATRIX_FLOAT(char* file, int start_row, int end_row, int total_column,
		MATRIX_FLOAT* mat, char* out_file);
int read_write_part_of_MATRIX_FLOAT(char* file, int start_row, int end_row, int total_column, char* out_file);
int read_write_part_of_MATRIX_FLOAT_alt(FILE* stream, int cnt_row, char* out_file);
void write_MATRIX_FLOAT(char* file, MATRIX_FLOAT mat);
void append_MATRIX_FLOAT_text(char* file, MATRIX_FLOAT mat);
void put_MATRIX_FLOAT(FILE* stream, MATRIX_FLOAT mat);
float norm1_float(float* v, int len);
void norm_rows_MATRIX_FLOAT(MATRIX_FLOAT mat, float p, VEC_FLOAT* norms);

void partitions( int n, int part_size, VEC_INT* part_sizes );
void partitions2( int n, int n_part, VEC_INT* part_sizes );
void calc_partitions_by_mem_size(int cnt_rows, int total_column, float max_mem_size,
		VEC_INT* part_sizes);
int quick_check_partitions_by_mem_size(int cnt_rows, int total_column, float max_mem_size, float *x);
#endif /* MATUTILS_H_ */
