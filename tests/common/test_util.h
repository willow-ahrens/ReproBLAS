#ifndef __TEST_VEC_H
#define __TEST_VEC_H

#include "../../src/types.h"
#include "../../build/include/rblas.h"

typedef enum util_vec_fill {
  util_Vec_Constant = 0,
  util_Vec_Rand,
  util_Vec_2_Times_Rand_Minus_1,
  util_Vec_Rand_Plus_Rand_Minus_1,
  util_Vec_Normal,
  util_Vec_Sine,
  util_Vec_Small_Plus_Increasing_Big,
  util_Vec_Small_Plus_Rand_Big,
  util_Vec_Rand_Cond,
  util_Vec_Constant_Drop,
  util_Vec_Rand_Drop,
  util_Vec_2_Times_Rand_Minus_1_Drop,
  util_Vec_Rand_Plus_Rand_Minus_1_Drop,
  util_Vec_Normal_Drop,
  util_Vec_Sine_Drop
} util_vec_fill_t;
static const int  util_vec_fill_n_names  = 15;
static const char *util_vec_fill_names[] = {"constant",
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
static const char *util_vec_fill_descs[] = {"Constant",
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
typedef enum util_mat_fill {
  util_Mat_Row_Constant = 0,
  util_Mat_Row_Rand,
  util_Mat_Row_2_Times_Rand_Minus_1,
  util_Mat_Row_Rand_Plus_Rand_Minus_1,
  util_Mat_Row_Normal,
  util_Mat_Row_Sine,
  util_Mat_Row_Small_Plus_Increasing_Big,
  util_Mat_Row_Small_Plus_Rand_Big,
  util_Mat_Row_Rand_Cond,
  util_Mat_Row_Constant_Drop,
  util_Mat_Row_Rand_Drop,
  util_Mat_Row_2_Times_Rand_Minus_1_Drop,
  util_Mat_Row_Rand_Plus_Rand_Minus_1_Drop,
  util_Mat_Row_Normal_Drop,
  util_Mat_Row_Sine_Drop,
  util_Mat_Identity
} util_mat_fill_t;
static const int  util_mat_fill_n_names  = 16;
static const char *util_mat_fill_names[] = {"constant",
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
                                            "sine[drop]",
                                            "identity"};
static const char *util_mat_fill_descs[] = {"Constant",
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
                                            "Sine(2pi*(i/n))[drop]",
                                            "Identity"};

void util_dvec_fill(int N, double* v, int inc, util_vec_fill_t fill, double a, double b);
void util_svec_fill(int N, float* v, int inc, util_vec_fill_t fill, float a, float b);
void util_zvec_fill(int N, double complex* v, int inc, util_vec_fill_t fill, double complex a, double complex b);
void util_cvec_fill(int N, float complex* v, int inc, util_vec_fill_t fill, float complex a, float complex b);

void util_dmat_fill(rblas_order_t order, rblas_transpose_t transA, int M, int N, double* A, int ldA, util_mat_fill_t fill, double a, double b);
void util_smat_fill(rblas_order_t order, rblas_transpose_t transA, int M, int N, float* A, int ldA, util_mat_fill_t fill, float a, float b);
void util_zmat_fill(rblas_order_t order, rblas_transpose_t transA, int M, int N, double complex* A, int ldA, util_mat_fill_t fill, double complex a, double complex b);
void util_cmat_fill(rblas_order_t order, rblas_transpose_t transA, int M, int N, float complex* A, int ldA, util_mat_fill_t fill, float complex a, float complex b);

void util_random_seed(void);

typedef enum util_comp {
  util_Increasing = 0,
  util_Increasing_Magnitude,
  util_Decreasing,
  util_Decreasing_Magnitude
} util_comp_t;
int util_dcompare(void *a, void *b, util_comp_t comp);
int util_scompare(void *a, void *b, util_comp_t comp);
int util_zcompare(void *a, void *b, util_comp_t comp);
int util_ccompare(void *a, void *b, util_comp_t comp);
void util_svec_sort(int N, float* V, int incV, int *P, int incP, util_comp_t comp);
void util_cvec_sort(int N, float complex* V, int incV, int *P, int incP, util_comp_t comp);
void util_dvec_sort(int N, double* V, int incV, int *P, int incP, util_comp_t comp);
void util_zvec_sort(int N, double complex* V, int incV, int *P, int incP, util_comp_t comp);

void util_svec_shuffle(int N, float* V, int incV, int *P, int incP);
void util_cvec_shuffle(int N, float complex* V, int incV, int *P, int incP);
void util_dvec_shuffle(int N, double* V, int incV, int *P, int incP);
void util_zvec_shuffle(int N, double complex* V, int incV, int *P, int incP);

void util_svec_reverse(int N, float* V, int incV, int *P, int incP);
void util_dvec_reverse(int N, double* V, int incV, int *P, int incP);
void util_cvec_reverse(int N, float complex* V, int incV, int *P, int incP);
void util_zvec_reverse(int N, double complex* V, int incV, int *P, int incP);

float* util_svec_alloc(int N, int incV);
double* util_dvec_alloc(int N, int incV);
float complex* util_cvec_alloc(int N, int incV);
double complex* util_zvec_alloc(int N, int incV);

double* util_dmat_alloc(rblas_order_t order, int M, int N, int ldA);
float* util_smat_alloc(rblas_order_t order, int M, int N, int ldA);
double complex* util_zmat_alloc(rblas_order_t order, int M, int N, int ldA);
float complex* util_cmat_alloc(rblas_order_t order, int M, int N, int ldA);


int* util_identity_permutation(int N);
int* util_inverse_permutation(int N, int *P, int incP);
void util_dvec_permute(int N, double* V, int incV, int *Q, int incQ, int *P, int incP);
void util_svec_permute(int N, float* V, int incV, int *Q, int incQ, int *P, int incP);
void util_zvec_permute(int N, double complex* V, int incV, int *Q, int incQ, int *P, int incP);
void util_cvec_permute(int N, float complex* V, int incV, int *Q, int incQ, int *P, int incP);

#endif
