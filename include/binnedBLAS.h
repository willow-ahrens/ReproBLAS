/**
 * @file  binnedBLAS.h
 * @brief binnedBLAS.h defines BLAS Methods that operate on binned types.
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
 * - db - binned double (#double_binned)
 * - zb - binned complex double (#double_complex_binned)
 * - sb - binned float (#float_binned)
 * - cb - binned complex float (#float_complex_binned)
 * - dm - manually specified binned double (@c double, @c double)
 * - zm - manually specified binned complex double (@c double, @c double)
 * - sm - manually specified binned float (@c float, @c float)
 * - cm - manually specified binned complex float (@c float, @c float)
 *
 * Throughout the library, complex types are specified via @c *void pointers. These routines will sometimes be suffixed by sub, to represent that a function has been made into a subroutine. This allows programmers to use whatever complex types they are already using, as long as the memory pointed to is of the form of two adjacent floating point types, the first and second representing real and imaginary components of the complex number.
 *
 * The goal of using binned types is to obtain either more accurate or reproducible summation of floating point numbers. In reproducible summation, floating point numbers are split into several slices along predefined boundaries in the exponent range. The space between two boundaries is called a bin. Indexed types are composed of several accumulators, each accumulating the slices in a particular bin. The accumulators correspond to the largest consecutive nonzero bins seen so far.
 *
 * The parameter @c fold describes how many accumulators are used in the binned types supplied to a subroutine (an binned type with @c k accumulators  is @c k-fold). The default value for this parameter can be set in config.h. If you are unsure of what value to use for @c fold, we recommend 3. Note that the @c fold of binned types must be the same for all binned types that interact with each other. Operations on more than one binned type assume all binned types being operated upon have the same @c fold. Note that the @c fold of an binned type may not be changed once the type has been allocated. A common use case would be to set the value of @c fold as a global macro in your code and supply it to all binned functions that you use.
 *
 * @internal
 * Power users of the library may find themselves wanting to manually specify the underlying primary and carry vectors of an binned type themselves. If you do not know what these are, don't worry about the manually specified binned types.
 */
#ifndef BINNEDBLAS_H_
#define BINNEDBLAS_H_
#include "binned.h"
#include "reproBLAS.h"

float binnedBLAS_samax(const int N, const float *X, const int incX);
double binnedBLAS_damax(const int N, const double *X, const int incX);
void binnedBLAS_camax_sub(const int N, const void *X, const int incX, void *amax);
void binnedBLAS_zamax_sub(const int N, const void *X, const int incX, void *amax);

float binnedBLAS_samaxm(const int N, const float *X, const int incX, const float *Y, const int incY);
double binnedBLAS_damaxm(const int N, const double *X, const int incX, const double *Y, const int incY);
void binnedBLAS_camaxm_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *amaxm);
void binnedBLAS_zamaxm_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *amaxm);

void binnedBLAS_dbdsum(const int fold, const int N, const double *X, const int incX, double_binned *Y);
void binnedBLAS_dmdsum(const int fold, const int N, const double *X, const int incX, double *priY, const int incpriY, double *carY, const int inccarY);
void binnedBLAS_dbdasum(const int fold, const int N, const double *X, const int incX, double_binned *Y);
void binnedBLAS_dmdasum(const int fold, const int N, const double *X, const int incX, double *priY, const int incpriY, double *carY, const int inccarY);
double binnedBLAS_dbdssq(const int fold, const int N, const double *X, const int incX, const double scaleY, double_binned *Y);
double binnedBLAS_dmdssq(const int fold, const int N, const double *X, const int incX, const double scaleY, double *priY, const int incpriY, double *carY, const int inccarY);
void binnedBLAS_dbddot(const int fold, const int N, const double *X, const int incX, const double *Y, const int incY, double_binned *Z);
void binnedBLAS_dmddot(const int fold, const int N, const double *X, const int incX, const double *Y, const int incY, double *manZ, const int incmanZ, double *carZ, const int inccarZ);

void binnedBLAS_zbzsum(const int fold, const int N, const void *X, const int incX, double_binned *Y);
void binnedBLAS_zmzsum(const int fold, const int N, const void *X, const int incX, double *priY, const int incpriY, double *carY, const int inccarY);
void binnedBLAS_dbzasum(const int fold, const int N, const void *X, const int incX, double_binned *Y);
void binnedBLAS_dmzasum(const int fold, const int N, const void *X, const int incX, double *priY, const int incpriY, double *carY, const int inccarY);
double binnedBLAS_dbzssq(const int fold, const int N, const void *X, const int incX, const double scaleY, double_binned *Y);
double binnedBLAS_dmzssq(const int fold, const int N, const void *X, const int incX, const double scaleY, double *priY, const int incpriY, double *carY, const int inccarY);
void binnedBLAS_zbzdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double_binned *Z);
void binnedBLAS_zmzdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double *manZ, const int incmanZ, double *carZ, const int inccarZ);
void binnedBLAS_zbzdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double_binned *Z);
void binnedBLAS_zmzdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double *manZ, const int incmanZ, double *carZ, const int inccarZ);

void binnedBLAS_sbssum(const int fold, const int N, const float *X, const int incX, float_binned *Y);
void binnedBLAS_smssum(const int fold, const int N, const float *X, const int incX, float *priY, const int incpriY, float *carY, const int inccarY);
void binnedBLAS_sbsasum(const int fold, const int N, const float *X, const int incX, float_binned *Y);
void binnedBLAS_smsasum(const int fold, const int N, const float *X, const int incX, float *priY, const int incpriY, float *carY, const int inccarY);
float binnedBLAS_sbsssq(const int fold, const int N, const float *X, const int incX, const float scaleY, float_binned *Y);
float binnedBLAS_smsssq(const int fold, const int N, const float *X, const int incX, const float scaleY, float *priY, const int incpriY, float *carY, const int inccarY);
void binnedBLAS_sbsdot(const int fold, const int N, const float *X, const int incX, const float *Y, const int incY, float_binned *Z);
void binnedBLAS_smsdot(const int fold, const int N, const float *X, const int incX, const float *Y, const int incY, float *manZ, const int incmanZ, float *carZ, const int inccarZ);

void binnedBLAS_cbcsum(const int fold, const int N, const void *X, const int incX, float_binned *Y);
void binnedBLAS_cmcsum(const int fold, const int N, const void *X, const int incX, float *priY, const int incpriY, float *carY, const int inccarY);
void binnedBLAS_sbcasum(const int fold, const int N, const void *X, const int incX, float_binned *Y);
void binnedBLAS_smcasum(const int fold, const int N, const void *X, const int incX, float *priY, const int incpriY, float *carY, const int inccarY);
float binnedBLAS_sbcssq(const int fold, const int N, const void *X, const int incX, const float scaleY, float_binned *Y);
float binnedBLAS_smcssq(const int fold, const int N, const void *X, const int incX, const float scaleY, float *priY, const int incpriY, float *carY, const int inccarY);
void binnedBLAS_cbcdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float_binned *Z);
void binnedBLAS_cmcdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float *manZ, const int incmanZ, float *carZ, const int inccarZ);
void binnedBLAS_cbcdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float_binned *Z);
void binnedBLAS_cmcdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float *manZ, const int incmanZ, float *carZ, const int inccarZ);

void binnedBLAS_dbdgemv(const int fold, const char Order, const char TransA,
             const int M, const int N,
             const double alpha, const double *A, const int lda,
             const double *X, const int incX,
             double_binned *Y, const int incY);
void binnedBLAS_dbdgemm(const int fold, const char Order,
             const char TransA, const char TransB,
             const int M, const int N, const int K,
             const double alpha, const double *A, const int lda,
             const double *B, const int ldb,
             double_binned *C, const int ldc);

void binnedBLAS_sbsgemv(const int fold, const char Order, const char TransA,
             const int M, const int N,
             const float alpha, const float *A, const int lda,
             const float *X, const int incX,
             float_binned *Y, const int incY);
void binnedBLAS_sbsgemm(const int fold, const char Order,
             const char TransA, const char TransB,
             const int M, const int N, const int K,
             const float alpha, const float *A, const int lda,
             const float *B, const int ldb,
             float_binned *C, const int ldc);

void binnedBLAS_zbzgemv(const int fold, const char Order, const char TransA,
             const int M, const int N,
             const void *alpha, const void *A, const int lda,
             const void *X, const int incX,
             double_complex_binned *Y, const int incY);
void binnedBLAS_zbzgemm(const int fold, const char Order,
             const char TransA, const char TransB,
             const int M, const int N, const int K,
             const void *alpha, const void *A, const int lda,
             const void *B, const int ldb,
             double_complex_binned *C, const int ldc);

void binnedBLAS_cbcgemv(const int fold, const char Order, const char TransA,
             const int M, const int N,
             const void *alpha, const void *A, const int lda,
             const void *X, const int incX,
             float_complex_binned *Y, const int incY);
void binnedBLAS_cbcgemm(const int fold, const char Order,
             const char TransA, const char TransB,
             const int M, const int N, const int K,
             const void *alpha, const void *A, const int lda,
             const void *B, const int ldb,
             float_complex_binned *C, const int ldc);

#endif
