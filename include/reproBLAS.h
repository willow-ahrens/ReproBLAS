/**
 * @file  reproBLAS.h
 * @brief reproBLAS.h defines reproducible BLAS Methods.
 *
 * This header is modeled after cblas.h, and as such functions are prefixed with character sets describing the data types they operate upon. For example, the function @c dfoo would perform the function @c foo on @c double possibly returning a @c double.
 *
 * If two character sets are prefixed, the first set of characters describes the output and the second the input type. For example, the function @c dzbar would perform the function @c bar on @c double @c complex and return a @c double.
 *
 * Such character sets are listed as follows:
 * - d - double (@c double)
 * - z - complex double (@c *void)
 * - s - float (@c float)
 * - c - complex float (@c *void)
 *
 * Throughout the library, complex types are specified via @c *void pointers. These routines will sometimes be suffixed by sub, to represent that a function has been made into a subroutine. This allows programmers to use whatever complex types they are already using, as long as the memory pointed to is of the form of two adjacent floating point types, the first and second representing real and imaginary components of the complex number.
 *
 * The goal of using binned types is to obtain either more accurate or reproducible summation of floating point numbers. In reproducible summation, floating point numbers are split into several slices along predefined boundaries in the exponent range. The space between two boundaries is called a bin. Binned types are composed of several accumulators, each accumulating the slices in a particular bin. The accumulators correspond to the largest consecutive nonzero bins seen so far.
 *
 * The parameter @c fold describes how many accumulators are used in the binned types supplied to a subroutine (an binned type with @c k accumulators  is @c k-fold). The default value for this parameter can be set in config.h. If you are unsure of what value to use for @c fold, we recommend 3. Note that the @c fold of binned types must be the same for all binned types that interact with each other. Operations on more than one binned type assume all binned types being operated upon have the same @c fold. Note that the @c fold of an binned type may not be changed once the type has been allocated. A common use case would be to set the value of @c fold as a global macro in your code and supply it to all binned functions that you use.
 *
 * In reproBLAS, two copies of the BLAS are provided. The functions that share the same name as their BLAS counterparts perform reproducible versions of their corresponding operations using the default fold value specified in config.h. The functions that are prefixed by the character 'r' allow the user to specify their own fold for the underlying binned types.
 */
#ifndef REPROBLAS_H_
#define REPROBLAS_H_
#include <complex.h>

double reproBLAS_rdsum(const int fold, const int N, const double* X, const int incX);
double reproBLAS_rdasum(const int fold, const int N, const double* X, const int incX);
double reproBLAS_rdnrm2(const int fold, const int N, const double* X, const int incX);
double reproBLAS_rddot(const int fold, const int N, const double* X, const int incX, const double* Y, const int incY);

float reproBLAS_rsdot(const int fold, const int N, const float* X, const int incX, const float* Y, const int incY);
float reproBLAS_rsasum(const int fold, const int N, const float* X, const int incX);
float reproBLAS_rssum(const int fold, const int N, const float* X, const int incX);
float reproBLAS_rsnrm2(const int fold, const int N, const float* X, const int incX);

void reproBLAS_rzsum_sub(const int fold, const int N, const void* X, int incX, void *sum);
double reproBLAS_rdzasum(const int fold, const int N, const void* X, const int incX);
double reproBLAS_rdznrm2(const int fold, const int N, const void* X, int incX);
void reproBLAS_rzdotc_sub(const int fold, const int N, const void* X, const int incX, const void* Y, const int incY, void *dotc);
void reproBLAS_rzdotu_sub(const int fold, const int N, const void* X, const int incX, const void* Y, const int incY, void *dotu);

void reproBLAS_rcsum_sub(const int fold, const int N, const void* X, const int incX, void *sum);
float reproBLAS_rscasum(const int fold, const int N, const void* X, const int incX);
float reproBLAS_rscnrm2(const int fold, const int N, const void* X, const int incX);
void reproBLAS_rcdotc_sub(const int fold, const int N, const void* X, const int incX, const void* Y, const int incY, void *dotc);
void reproBLAS_rcdotu_sub(const int fold, const int N, const void* X, const int incX, const void* Y, const int incY, void *dotu);

void reproBLAS_rdgemv(const int fold, const char Order, const char TransA,
            const int M, const int N,
            const double alpha, const double *A, const int lda,
            const double *X, const int incX,
            const double beta, double *Y, const int incY);
void reproBLAS_rdgemm(const int fold, const char Order,
            const char TransA, const char TransB,
            const int M, const int N, const int K,
            const double alpha, const double *A, const int lda,
            const double *B, const int ldb,
            const double beta, double *C, const int ldc);

void reproBLAS_rsgemv(const int fold, const char Order, const char TransA,
            const int M, const int N,
            const float alpha, const float *A, const int lda,
            const float *X, const int incX,
            const float beta, float *Y, const int incY);
void reproBLAS_rsgemm(const int fold, const char Order,
            const char TransA, const char TransB,
            const int M, const int N, const int K,
            const float alpha, const float *A, const int lda,
            const float *B, const int ldb,
            const float beta, float *C, const int ldc);

void reproBLAS_rzgemv(const int fold, const char Order, const char TransA,
            const int M, const int N,
            const void *alpha, const void *A, const int lda,
            const void *X, const int incX,
            const void *beta, void *Y, const int incY);
void reproBLAS_rzgemm(const int fold, const char Order,
            const char TransA, const char TransB,
            const int M, const int N, const int K,
            const void *alpha, const void *A, const int lda,
            const void *B, const int ldb,
            const void *beta, void *C, const int ldc);

void reproBLAS_rcgemv(const int fold, const char Order, const char TransA,
            const int M, const int N,
            const void *alpha, const void *A, const int lda,
            const void *X, const int incX,
            const void *beta, void *Y, const int incY);
void reproBLAS_rcgemm(const int fold, const char Order,
            const char TransA, const char TransB,
            const int M, const int N, const int K,
            const void *alpha, const void *A, const int lda,
            const void *B, const int ldb,
            const void *beta, void *C, const int ldc);

double reproBLAS_dsum(const int N, const double* X, const int incX);
double reproBLAS_dasum(const int N, const double* X, const int incX);
double reproBLAS_dnrm2(const int N, const double* X, const int incX);
double reproBLAS_ddot(const int N, const double* X, const int incX, const double* Y, const int incY);

float reproBLAS_sdot(const int N, const float* X, const int incX, const float* Y, const int incY);
float reproBLAS_sasum(const int N, const float* X, const int incX);
float reproBLAS_ssum(const int N, const float* X, const int incX);
float reproBLAS_snrm2(const int N, const float* X, const int incX);

void reproBLAS_zsum_sub(const int N, const void* X, int incX, void *sum);
double reproBLAS_dzasum(const int N, const void* X, const int incX);
double reproBLAS_dznrm2(const int N, const void* X, int incX);
void reproBLAS_zdotc_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotc);
void reproBLAS_zdotu_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotu);

void reproBLAS_csum_sub(const int N, const void* X, const int incX, void *sum);
float reproBLAS_scasum(const int N, const void* X, const int incX);
float reproBLAS_scnrm2(const int N, const void* X, const int incX);
void reproBLAS_cdotc_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotc);
void reproBLAS_cdotu_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotu);

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

void reproBLAS_sgemv(const char Order, const char TransA,
            const int M, const int N,
            const float alpha, const float *A, const int lda,
            const float *X, const int incX,
            const float beta, float *Y, const int incY);
void reproBLAS_sgemm(const char Order, const char TransA, const char TransB,
            const int M, const int N, const int K,
            const float alpha, const float *A, const int lda,
            const float *B, const int ldb,
            const float beta, float *C, const int ldc);

void reproBLAS_zgemv(const char Order, const char TransA,
            const int M, const int N,
            const void *alpha, const void *A, const int lda,
            const void *X, const int incX,
            const void *beta, void *Y, const int incY);
void reproBLAS_zgemm(const char Order, const char TransA, const char TransB,
            const int M, const int N, const int K,
            const void *alpha, const void *A, const int lda,
            const void *B, const int ldb,
            const void *beta, void *C, const int ldc);

void reproBLAS_cgemv(const char Order, const char TransA,
            const int M, const int N,
            const void *alpha, const void *A, const int lda,
            const void *X, const int incX,
            const void *beta, void *Y, const int incY);
void reproBLAS_cgemm(const char Order, const char TransA, const char TransB,
            const int M, const int N, const int K,
            const void *alpha, const void *A, const int lda,
            const void *B, const int ldb,
            const void *beta, void *C, const int ldc);

#endif
