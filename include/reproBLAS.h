#ifndef REPROBLAS_H_
#define REPROBLAS_H_
#include <complex.h>

double reproBLAS_rdsum(const int fold, const int N, const double* X, const int incX);
double reproBLAS_rdasum(const int fold, const int N, const double* X, const int incX);
double reproBLAS_rdnrm2(const int fold, const int N, const double* X, const int incX);
double reproBLAS_rddot(const int fold, const int N, const double* X, const int incX, const double* Y, const int incY);
void reproBLAS_rdgemv(const int fold, const char Order, const char TransA,
            const int M, const int N,
            const double alpha, const double *A, const int lda,
            const double *X, const int incX,
            const double beta, double *Y, const int incY);
void reproBLAS_rdgemm(const int fold, const char Order, const char TransA, const char TransB,
            const int M, const int N, const int K,
            const double alpha, const double *A, const int lda,
            const double *B, const int ldb,
            const double beta, double *C, const int ldc);

float reproBLAS_rsdot(const int fold, const int N, const float* X, const int incX, const float* Y, const int incY);
float reproBLAS_rsasum(const int fold, const int N, const float* X, const int incX);
float reproBLAS_rssum(const int fold, const int N, const float* X, const int incX);
float reproBLAS_rsnrm2(const int fold, const int N, const float* X, const int incX);

void reproBLAS_rzsum_sub(const int fold, const int N, const void* X, int incX, void *sum);
double reproBLAS_rdzasum(const int fold, const int N, const void* X, const int incX);
double reproBLAS_rdznrm2(const int fold, const int N, const void* X, int incX);
void reproBLAS_rzdotc_sub(const int fold, const int N, const void* X, const int incX, const void* Y, const int incY, void *dotc);
void reproBLAS_rzdotu_sub(const int fold, const int N, const void* X, const int incX, const void* Y, const int incY, void *dotu);
void reproBLAS_rzgemv(const int fold, const char Order, const char TransA,
            const int M, const int N,
            const void *alpha, const void *A, const int lda,
            const void *X, const int incX,
            const void *beta, void *Y, const int incY);
void reproBLAS_rzgemm(const int fold, const char Order, const char TransA, const char TransB,
            const int M, const int N, const int K,
            const void *alpha, const void *A, const int lda,
            const void *B, const int ldb,
            const void *beta, void *C, const int ldc);

void reproBLAS_rcsum_sub(const int fold, const int N, const void* X, const int incX, void *sum);
float reproBLAS_rscasum(const int fold, const int N, const void* X, const int incX);
float reproBLAS_rscnrm2(const int fold, const int N, const void* X, const int incX);
void reproBLAS_rcdotc_sub(const int fold, const int N, const void* X, const int incX, const void* Y, const int incY, void *dotc);
void reproBLAS_rcdotu_sub(const int fold, const int N, const void* X, const int incX, const void* Y, const int incY, void *dotu);

double reproBLAS_dsum(const int N, const double* X, const int incX);
double reproBLAS_dasum(const int N, const double* X, const int incX);
double reproBLAS_dnrm2(const int N, const double* X, const int incX);
double reproBLAS_ddot(const int N, const double* X, const int incX, const double* Y, const int incY);
void reproBLAS_dgemv(const char Order, const char TransA,
            const int M, const int N,
            const double alpha, const double *A, const int lda,
            const double *X, const int incX,
            const double beta, double *Y, const int incY);
void reproBLAS_dgemm(const char Order, const char TransA, const char TransB,
            const int M, const int N, const int K,
            const double alpha, const double *A, const int lda,
            const double *B, const int ldb,
            const double beta, double *C, const int ldc);

float reproBLAS_sdot(const int N, const float* X, const int incX, const float* Y, const int incY);
float reproBLAS_sasum(const int N, const float* X, const int incX);
float reproBLAS_ssum(const int N, const float* X, const int incX);
float reproBLAS_snrm2(const int N, const float* X, const int incX);

void reproBLAS_zsum_sub(const int N, const void* X, int incX, void *sum);
double reproBLAS_dzasum(const int N, const void* X, const int incX);
double reproBLAS_dznrm2(const int N, const void* X, int incX);
void reproBLAS_zdotc_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotc);
void reproBLAS_zdotu_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotu);
void reproBLAS_zgemv(const char Order, const char TransA,
            const int M, const int N,
            const void *alpha, const void *A, const int lda,
            const void *X, const int incX,
            const void *beta, void *Y, const int incY);

void reproBLAS_csum_sub(const int N, const void* X, const int incX, void *sum);
float reproBLAS_scasum(const int N, const void* X, const int incX);
float reproBLAS_scnrm2(const int N, const void* X, const int incX);
void reproBLAS_cdotc_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotc);
void reproBLAS_cdotu_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotu);

#endif
