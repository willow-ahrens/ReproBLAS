#ifndef TEST_UTIL_H
#define TEST_UTIL_H

#include <complex.h>

typedef enum util_vec_fill {
  util_Vec_Constant = 0,
  util_Vec_Pos_Inf,
  util_Vec_Pos_Pos_Inf,
  util_Vec_Pos_Neg_Inf,
  util_Vec_NaN,
  util_Vec_Pos_Inf_NaN,
  util_Vec_Pos_Pos_Inf_NaN,
  util_Vec_Pos_Neg_Inf_NaN,
  util_Vec_Pos_Big,
  util_Vec_Pos_Pos_Big,
  util_Vec_Pos_Neg_Big,
  util_Vec_Rand,
  util_Vec_2_Times_Rand_Minus_1,
  util_Vec_Rand_Plus_Rand_Minus_1,
  util_Vec_Normal,
  util_Vec_Sine,
  util_Vec_Small_Plus_Increasing_Big,
  util_Vec_Small_Plus_Rand_Big,
  util_Vec_Rand_Cond,
  util_Vec_N_Cond,
  util_Vec_N_Cond_Plus_Rand,
  util_Vec_Constant_Drop,
  util_Vec_Rand_Drop,
  util_Vec_2_Times_Rand_Minus_1_Drop,
  util_Vec_Rand_Plus_Rand_Minus_1_Drop,
  util_Vec_Normal_Drop,
  util_Vec_Sine_Drop
} util_vec_fill_t;
extern const int util_vec_fill_n_names;
extern const char *util_vec_fill_names[];
extern const char *util_vec_fill_descs[];

typedef enum util_mat_fill {
  util_Mat_Row_Constant = 0,
  util_Mat_Row_Pos_Inf,
  util_Mat_Row_Pos_Pos_Inf,
  util_Mat_Row_Pos_Neg_Inf,
  util_Mat_Row_NaN,
  util_Mat_Row_Pos_Inf_NaN,
  util_Mat_Row_Pos_Pos_Inf_NaN,
  util_Mat_Row_Pos_Neg_Inf_NaN,
  util_Mat_Row_Pos_Big,
  util_Mat_Row_Pos_Pos_Big,
  util_Mat_Row_Pos_Neg_Big,
  util_Mat_Row_Rand,
  util_Mat_Row_2_Times_Rand_Minus_1,
  util_Mat_Row_Rand_Plus_Rand_Minus_1,
  util_Mat_Row_Normal,
  util_Mat_Row_Sine,
  util_Mat_Row_Small_Plus_Increasing_Big,
  util_Mat_Row_Small_Plus_Rand_Big,
  util_Mat_Row_Rand_Cond,
  util_Mat_Row_N_Cond,
  util_Mat_Row_N_Cond_Plus_Rand,
  util_Mat_Row_Constant_Drop,
  util_Mat_Row_Rand_Drop,
  util_Mat_Row_2_Times_Rand_Minus_1_Drop,
  util_Mat_Row_Rand_Plus_Rand_Minus_1_Drop,
  util_Mat_Row_Normal_Drop,
  util_Mat_Row_Sine_Drop,
  util_Mat_Identity
} util_mat_fill_t;
extern const int  util_mat_fill_n_names;
extern const char *util_mat_fill_names[];
extern const char *util_mat_fill_descs[];

double util_drand48();

void util_ddpd(double* a, double b);

#ifndef MIN
  #define MIN(A, B) (((A) < (B))? (A): (B))
#endif

#ifndef MAX
  #define MAX(A, B) (((A) > (B))? (A): (B))
#endif

int util_dsoftequals(double a, double b, double bound);
int util_zsoftequals(double complex a, double complex b, double complex bound);
int util_ssoftequals(float a, float b, float bound);
int util_csoftequals(float complex a, float complex b, float complex bound);

void util_dvec_fill(int N, double* v, int inc, util_vec_fill_t Fill, double RealScale, double ImagScale);
void util_svec_fill(int N, float* v, int inc, util_vec_fill_t Fill, float RealScale, float ImagScale);
void util_zvec_fill(int N, double complex* v, int inc, util_vec_fill_t Fill, double RealScale, double ImagScale);
void util_cvec_fill(int N, float complex* v, int inc, util_vec_fill_t Fill, float RealScale, float ImagScale);

void util_dmat_fill(char Order, char TransA, int M, int N, double* A, int lda, util_mat_fill_t Fill, double RealScale, double ImagScale);
void util_smat_fill(char Order, char TransA, int M, int N, float* A, int lda, util_mat_fill_t Fill, float RealScale, float ImagScale);
void util_zmat_fill(char Order, char TransA, int M, int N, double complex* A, int lda, util_mat_fill_t Fill, double RealScale, double ImagScale);
void util_cmat_fill(char Order, char TransA, int M, int N, float complex* A, int lda, util_mat_fill_t Fill, float RealScale, float ImagScale);

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
void util_dmat_row_sort(char Order, char TransA, int M, int N, double *A, int lda, int *P, int incP, util_comp_t comp, int col);
void util_smat_row_sort(char Order, char TransA, int M, int N, float *A, int lda, int *P, int incP, util_comp_t comp, int col);
void util_zmat_row_sort(char Order, char TransA, int M, int N, double complex *A, int lda, int *P, int incP, util_comp_t comp, int col);
void util_cmat_row_sort(char Order, char TransA, int M, int N, float complex *A, int lda, int *P, int incP, util_comp_t comp, int col);

void util_svec_shuffle(int N, float* V, int incV, int *P, int incP);
void util_cvec_shuffle(int N, float complex* V, int incV, int *P, int incP);
void util_dvec_shuffle(int N, double* V, int incV, int *P, int incP);
void util_zvec_shuffle(int N, double complex* V, int incV, int *P, int incP);
void util_dmat_row_shuffle(char Order, char TransA, int M, int N, double *A, int lda, int *P, int incP);
void util_smat_row_shuffle(char Order, char TransA, int M, int N, float *A, int lda, int *P, int incP);
void util_zmat_row_shuffle(char Order, char TransA, int M, int N, double complex *A, int lda, int *P, int incP);
void util_cmat_row_shuffle(char Order, char TransA, int M, int N, float complex *A, int lda, int *P, int incP);

void util_svec_reverse(int N, float* V, int incV, int *P, int incP);
void util_dvec_reverse(int N, double* V, int incV, int *P, int incP);
void util_cvec_reverse(int N, float complex* V, int incV, int *P, int incP);
void util_zvec_reverse(int N, double complex* V, int incV, int *P, int incP);
void util_dmat_row_reverse(char Order, char TransA, int M, int N, double *A, int lda, int *P, int incP);
void util_smat_row_reverse(char Order, char TransA, int M, int N, float *A, int lda, int *P, int incP);
void util_zmat_row_reverse(char Order, char TransA, int M, int N, double complex *A, int lda, int *P, int incP);
void util_cmat_row_reverse(char Order, char TransA, int M, int N, float complex *A, int lda, int *P, int incP);

float* util_svec_alloc(int N, int incV);
double* util_dvec_alloc(int N, int incV);
float complex* util_cvec_alloc(int N, int incV);
double complex* util_zvec_alloc(int N, int incV);

double* util_dmat_alloc(char Order, int M, int N, int lda);
float* util_smat_alloc(char Order, int M, int N, int lda);
double complex* util_zmat_alloc(char Order, int M, int N, int lda);
float complex* util_cmat_alloc(char Order, int M, int N, int lda);

int* util_identity_permutation(int N);
int* util_inverse_permutation(int N, int *P, int incP);
void util_dvec_permute(int N, double* V, int incV, int *Q, int incQ, int *P, int incP);
void util_svec_permute(int N, float* V, int incV, int *Q, int incQ, int *P, int incP);
void util_zvec_permute(int N, double complex* V, int incV, int *Q, int incQ, int *P, int incP);
void util_cvec_permute(int N, float complex* V, int incV, int *Q, int incQ, int *P, int incP);
void util_divec_permute(int num, int N, double* V, int incV, int *Q, int incQ, int *P, int incP);
void util_sivec_permute(int num, int N, float* V, int incV, int *Q, int incQ, int *P, int incP);
void util_zivec_permute(int num, int N, double complex* V, int incV, int *Q, int incQ, int *P, int incP);
void util_civec_permute(int num, int N, float complex* V, int incV, int *Q, int incQ, int *P, int incP);
void util_dmat_row_permute(char Order, char TransA, int M, int N, double *A, int lda, int *Q, int incQ, int *P, int incP);
void util_smat_row_permute(char Order, char TransA, int M, int N, float *A, int lda, int *Q, int incQ, int *P, int incP);
void util_zmat_row_permute(char Order, char TransA, int M, int N, double complex *A, int lda, int *Q, int incQ, int *P, int incP);
void util_cmat_row_permute(char Order, char TransA, int M, int N, float complex *A, int lda, int *Q, int incQ, int *P, int incP);
void util_dimat_row_permute(int num, char Order, char TransA, int M, int N, double *A, int lda, int *Q, int incQ, int *P, int incP);
void util_simat_row_permute(int num, char Order, char TransA, int M, int N, float *A, int lda, int *Q, int incQ, int *P, int incP);
void util_zimat_row_permute(int num, char Order, char TransA, int M, int N, double complex *A, int lda, int *Q, int incQ, int *P, int incP);
void util_cimat_row_permute(int num, char Order, char TransA, int M, int N, float complex *A, int lda, int *Q, int incQ, int *P, int incP);

#endif
