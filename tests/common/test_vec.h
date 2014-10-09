#include "../../src/types.h"

#define vec_fill_CONSTANT                    0
#define vec_fill_RAND                        1
#define vec_fill_2_TIMES_RAND_MINUS_1        2
#define vec_fill_RAND_PLUS_RAND_MINUS_1      3
#define vec_fill_NORMAL                      4
#define vec_fill_SINE                        5
#define vec_fill_SMALL_PLUS_INCREASING_BIG   6
#define vec_fill_SMALL_PLUS_RAND_BIG         7
#define vec_fill_RAND_COND                   8
#define vec_fill_CONSTANT_DROP               9
#define vec_fill_RAND_DROP                   10
#define vec_fill_2_TIMES_RAND_MINUS_1_DROP   11
#define vec_fill_RAND_PLUS_RAND_MINUS_1_DROP 12
#define vec_fill_NORMAL_DROP                 13
#define vec_fill_SINE_DROP                   14
#define vec_fill_MAX                         15
void svec_fill(int N, float* v, int inc, int type, float a, float b);
void dvec_fill(int N, double* v, int inc, int type, double a, double b);
void cvec_fill(int N, float complex* v, int inc, int type, float complex a, float complex b);
void zvec_fill(int N, double complex* v, int inc, int type, double complex a, double complex b);

const char* vec_fill_name(int type);

void vec_random_seed(void);

#define vec_order_INCREASING            0
#define vec_order_INCREASING_MAGNITUDE  1
#define vec_order_DECREASING            2
#define vec_order_DECREASING_MAGNITUDE  3
void svec_sort(int N, float* v, int inc, int type);
void cvec_sort(int N, float complex* v, int inc, int type);
void dvec_sort(int N, double* v, int inc, int type);
void zvec_sort(int N, double complex* v, int inc, int type);

void svec_shuffle(int N, float* v, int inc);
void cvec_shuffle(int N, float complex* v, int inc);
void dvec_shuffle(int N, double* v, int inc);
void zvec_shuffle(int N, double complex* v, int inc);

void svec_reverse(int N, float* v, int inc);
void dvec_reverse(int N, double* v, int inc);
void cvec_reverse(int N, float complex* v, int inc);
void zvec_reverse(int N, double complex* v, int inc);

float* svec_alloc(int N, int inc);
double* dvec_alloc(int N, int inc);
float complex* cvec_alloc(int N, int inc);
double complex* zvec_alloc(int N, int inc);
