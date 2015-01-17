#ifndef __TEST_VEC_H
#define __TEST_VEC_H

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
static const int  vec_fill_n_names  = 15;
static const char *vec_fill_names[] = {"constant",
                                       "rand",
                                       "2*rand-1",
                                       "rand+(rand-1)",
                                       "normal",
                                       "sine",
                                       "small+grow*big",
                                       "small+rand*big",
                                       "rand_cond",
                                       "constant[drop]",
                                       "rand[drop]",
                                       "2*rand-1[drop]",
                                       "rand+(rand-1)[drop]",
                                       "normal[drop]",
                                       "sine[drop]"};
static const char *vec_fill_descs[] = {"Constant",
                                       "Random",
                                       "2*Random-1",
                                       "Random+(Random-1)",
                                       "Normal",
                                       "Sine(2pi*(i/n))",
                                       "Small+(i/n)*Big",
                                       "Small+Rand*Big",
                                       "RandomConditioned",
                                       "Constant[drop]",
                                       "Random[drop]",
                                       "2*Random-1[drop]",
                                       "Random+(Random-1)[drop]",
                                       "Normal[drop]",
                                       "Sine(2pi*(i/n))[drop]"};

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
int dutil_compare(void *a, void *b, int order);
int sutil_compare(void *a, void *b, int order);
int zutil_compare(void *a, void *b, int order);
int cutil_compare(void *a, void *b, int order);
void svec_sort(int N, float* V, int incV, int *P, int incP, int order);
void cvec_sort(int N, float complex* V, int incV, int *P, int incP, int order);
void dvec_sort(int N, double* V, int incV, int *P, int incP, int order);
void zvec_sort(int N, double complex* V, int incV, int *P, int incP, int order);

void svec_shuffle(int N, float* V, int incV, int *P, int incP);
void cvec_shuffle(int N, float complex* V, int incV, int *P, int incP);
void dvec_shuffle(int N, double* V, int incV, int *P, int incP);
void zvec_shuffle(int N, double complex* V, int incV, int *P, int incP);

void svec_reverse(int N, float* V, int incV, int *P, int incP);
void dvec_reverse(int N, double* V, int incV, int *P, int incP);
void cvec_reverse(int N, float complex* V, int incV, int *P, int incP);
void zvec_reverse(int N, double complex* V, int incV, int *P, int incP);

float* svec_alloc(int N, int incV);
double* dvec_alloc(int N, int incV);
float complex* cvec_alloc(int N, int incV);
double complex* zvec_alloc(int N, int incV);

int* util_identity_permutation(int N);
int* util_inverse_permutation(int N, int *P, int incP);
void dvec_permute(int N, double* V, int incV, int *Q, int incQ, int *P, int incP);
void svec_permute(int N, float* V, int incV, int *Q, int incQ, int *P, int incP);
void zvec_permute(int N, double complex* V, int incV, int *Q, int incQ, int *P, int incP);
void cvec_permute(int N, float complex* V, int incV, int *Q, int incQ, int *P, int incP);

#endif
