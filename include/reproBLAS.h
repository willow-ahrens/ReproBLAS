#ifndef _REPRODUCIBLE_BLAS__H_
#define _REPRODUCIBLE_BLAS__H_
#include <complex.h>

double rdsum(const int N, const double* X, const int incX);
double rdasum(const int N, const double* X, const int incX);
double rdnrm2(const int N, const double* X, const int incX);
double rddot(const int N, const double* X, const int incX, const double* Y, const int incY);

float rsdot(const int N, const float* X, const int incX, const float* Y, const int incY);
float rsasum(const int N, const float* X, const int incX);
float rssum(const int N, const float* X, const int incX);
float rsnrm2(const int N, const float* X, const int incX);

void rzsum_sub(const int N, const void* X, int incX, void *sum);
double rdzasum(const int N, const void* X, const int incX);
double rdznrm2(const int N, const void* X, int incX);
void rzdotc_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotc);
void rzdotu_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotu);

void rcsum_sub(const int N, const void* X, const int incX, void *sum);
float rscasum(const int N, const void* X, const int incX);
float rscnrm2(const int N, const void* X, const int incX);
void rcdotc_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotc);
void rcdotu_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotu);

void rdgemv(const char Order, const char TransA,
            const int M, const int N,
            const double alpha, const double *A, const int lda,
            const double *X, const int incX,
            const double beta, double *Y, const int incY);
#endif
