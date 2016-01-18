/**
 * @file  idxdBLAS.h
 * @brief idxdBLAS.h defines BLAS Methods that operate on indexed types.
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
 * - di - indexed double (#double_indexed)
 * - zi - indexed complex double (#double_complex_indexed)
 * - si - indexed float (#float_indexed)
 * - ci - indexed complex float (#float_complex_indexed)
 * - dm - manually specified indexed double (@c double, @c double)
 * - zm - manually specified indexed complex double (@c double, @c double)
 * - sm - manually specified indexed float (@c float, @c float)
 * - cm - manually specified indexed complex float (@c float, @c float)
 *
 * Throughout the library, complex types are specified via @c *void pointers. These routines will sometimes be suffixed by sub, to represent that a function has been made into a subroutine. This allows programmers to use whatever complex types they are already using, as long as the memory pointed to is of the form of two adjacent floating point types, the first and second representing real and imaginary components of the complex number.
 *
 * The goal of using indexed types is to obtain either more accurate or reproducible summation of floating point numbers. Indexed types are composed of several adjacent bins...
 *
 * The parameter @c fold describes how many bins are used in the indexed types supplied to a subroutine. The maximum value for this parameter can be set in config.h. If you are unsure of what value to use for @fold, we recommend 3. Note that the @c fold of indexed types must be the same for all indexed types that interact with each other. Operations on more than one indexed type assume all indexed types being operated upon have the same @c fold. Note that the @c fold of an indexed type may not be changed once the type has been allocated. A common use case would be to set the value of @c fold as a global macro in your code and supply it to all indexed functions that you use.
 *
 * @internal
 * Power users of the library may find themselves wanting to manually specify the underlying primary and carry vectors of an indexed type themselves. If you do not know what these are, don't worry about the manually specified indexed types.
 */
#ifndef IDXDBLAS_H_
#define IDXDBLAS_H_
#include "idxd.h"
#include "reproBLAS.h"
#include <complex.h>

float idxdBLAS_samax(const int N, const float *X, const int incX);
double idxdBLAS_damax(const int N, const double *X, const int incX);
void idxdBLAS_camax_sub(const int N, const void *X, const int incX, void *amax);
void idxdBLAS_zamax_sub(const int N, const void *X, const int incX, void *amax);

float idxdBLAS_samaxm(const int N, const float *X, const int incX, const float *Y, const int incY);
double idxdBLAS_damaxm(const int N, const double *X, const int incX, const double *Y, const int incY);
void idxdBLAS_camaxm_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *amaxm);
void idxdBLAS_zamaxm_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *amaxm);

void idxdBLAS_didsum(const int fold, const int N, const double *X, const int incX, double_indexed *Y);
void idxdBLAS_dmdsum(const int fold, const int N, const double *X, const int incX, double *priY, const int incpriY, double *carY, const int inccarY);
void idxdBLAS_didasum(const int fold, const int N, const double *X, const int incX, double_indexed *Y);
void idxdBLAS_dmdasum(const int fold, const int N, const double *X, const int incX, double *priY, const int incpriY, double *carY, const int inccarY);
double idxdBLAS_didssq(const int fold, const int N, const double *X, const int incX, const double scaleY, double_indexed *Y);
double idxdBLAS_dmdssq(const int fold, const int N, const double *X, const int incX, const double scaleY, double *priY, const int incpriY, double *carY, const int inccarY);
void idxdBLAS_diddot(const int fold, const int N, const double *X, const int incX, const double *Y, const int incY, double_indexed *Z);
void idxdBLAS_dmddot(const int fold, const int N, const double *X, const int incX, const double *Y, const int incY, double *manZ, const int incmanZ, double *carZ, const int inccarZ);

void idxdBLAS_zizsum(const int fold, const int N, const void *X, const int incX, double_indexed *Y);
void idxdBLAS_zmzsum(const int fold, const int N, const void *X, const int incX, double *priY, const int incpriY, double *carY, const int inccarY);
void idxdBLAS_dizasum(const int fold, const int N, const void *X, const int incX, double_indexed *Y);
void idxdBLAS_dmzasum(const int fold, const int N, const void *X, const int incX, double *priY, const int incpriY, double *carY, const int inccarY);
double idxdBLAS_dizssq(const int fold, const int N, const void *X, const int incX, const double scaleY, double_indexed *Y);
double idxdBLAS_dmzssq(const int fold, const int N, const void *X, const int incX, const double scaleY, double *priY, const int incpriY, double *carY, const int inccarY);
void idxdBLAS_zizdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double_indexed *Z);
void idxdBLAS_zmzdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double *manZ, const int incmanZ, double *carZ, const int inccarZ);
void idxdBLAS_zizdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double_indexed *Z);
void idxdBLAS_zmzdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double *manZ, const int incmanZ, double *carZ, const int inccarZ);

void idxdBLAS_sissum(const int fold, const int N, const float *X, const int incX, float_indexed *Y);
void idxdBLAS_smssum(const int fold, const int N, const float *X, const int incX, float *priY, const int incpriY, float *carY, const int inccarY);
void idxdBLAS_sisasum(const int fold, const int N, const float *X, const int incX, float_indexed *Y);
void idxdBLAS_smsasum(const int fold, const int N, const float *X, const int incX, float *priY, const int incpriY, float *carY, const int inccarY);
float idxdBLAS_sisssq(const int fold, const int N, const float *X, const int incX, const float scaleY, float_indexed *Y);
float idxdBLAS_smsssq(const int fold, const int N, const float *X, const int incX, const float scaleY, float *priY, const int incpriY, float *carY, const int inccarY);
void idxdBLAS_sisdot(const int fold, const int N, const float *X, const int incX, const float *Y, const int incY, float_indexed *Z);
void idxdBLAS_smsdot(const int fold, const int N, const float *X, const int incX, const float *Y, const int incY, float *manZ, const int incmanZ, float *carZ, const int inccarZ);

void idxdBLAS_cicsum(const int fold, const int N, const void *X, const int incX, float_indexed *Y);
void idxdBLAS_cmcsum(const int fold, const int N, const void *X, const int incX, float *priY, const int incpriY, float *carY, const int inccarY);
void idxdBLAS_sicasum(const int fold, const int N, const void *X, const int incX, float_indexed *Y);
void idxdBLAS_smcasum(const int fold, const int N, const void *X, const int incX, float *priY, const int incpriY, float *carY, const int inccarY);
float idxdBLAS_sicssq(const int fold, const int N, const void *X, const int incX, const float scaleY, float_indexed *Y);
float idxdBLAS_smcssq(const int fold, const int N, const void *X, const int incX, const float scaleY, float *priY, const int incpriY, float *carY, const int inccarY);
void idxdBLAS_cicdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float_indexed *Z);
void idxdBLAS_cmcdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float *manZ, const int incmanZ, float *carZ, const int inccarZ);
void idxdBLAS_cicdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float_indexed *Z);
void idxdBLAS_cmcdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float *manZ, const int incmanZ, float *carZ, const int inccarZ);

void idxdBLAS_didgemv(const int fold, const char Order, const char TransA,
             const int M, const int N,
             const double alpha, const double *A, const int lda,
             const double *X, const int incX,
             double_indexed *Y, const int incY);
void idxdBLAS_didgemm(const int fold, const char Order,
             const char TransA, const char TransB,
             const int M, const int N, const int K,
             const double alpha, const double *A, const int lda,
             const double *B, const int ldb,
             double_indexed *C, const int ldc);

void idxdBLAS_sisgemv(const int fold, const char Order, const char TransA,
             const int M, const int N,
             const float alpha, const float *A, const int lda,
             const float *X, const int incX,
             float_indexed *Y, const int incY);
void idxdBLAS_sisgemm(const int fold, const char Order,
             const char TransA, const char TransB,
             const int M, const int N, const int K,
             const float alpha, const float *A, const int lda,
             const float *B, const int ldb,
             float_indexed *C, const int ldc);

void idxdBLAS_zizgemv(const int fold, const char Order, const char TransA,
             const int M, const int N,
             const void *alpha, const void *A, const int lda,
             const void *X, const int incX,
             double_complex_indexed *Y, const int incY);
void idxdBLAS_zizgemm(const int fold, const char Order,
             const char TransA, const char TransB,
             const int M, const int N, const int K,
             const void *alpha, const void *A, const int lda,
             const void *B, const int ldb,
             double_complex_indexed *C, const int ldc);

void idxdBLAS_cicgemv(const int fold, const char Order, const char TransA,
             const int M, const int N,
             const void *alpha, const void *A, const int lda,
             const void *X, const int incX,
             float_complex_indexed *Y, const int incY);
void idxdBLAS_cicgemm(const int fold, const char Order,
             const char TransA, const char TransB,
             const int M, const int N, const int K,
             const void *alpha, const void *A, const int lda,
             const void *B, const int ldb,
             float_complex_indexed *C, const int ldc);

#endif
