/**
 * @file  indexedBLAS.h
 * @brief indexedBLAS.h defines BLAS Methods that operate on indexed types.
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
#ifndef _INDEXED_BLAS__H_
#define _INDEXED_BLAS__H_
#include "indexed.h"
#include "reproBLAS.h"
#include <complex.h>

float samax(const int N, const float *X, const int incX);
double damax(const int N, const double *X, const int incX);
void camax_sub(const int N, const void *X, const int incX, void *amax);
void zamax_sub(const int N, const void *X, const int incX, void *amax);

float samaxm(const int N, const float *X, const int incX, const float *Y, const int incY);
double damaxm(const int N, const double *X, const int incX, const double *Y, const int incY);
void camaxm_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *amaxm);
void zamaxm_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *amaxm);

void didsum(const int fold, const int N, const double *X, const int incX, double_indexed *Y);
void dmdsum(const int fold, const int N, const double *X, const int incX, double *priY, const int incpriY, double *carY, const int inccarY);
void didasum(const int fold, const int N, const double *X, const int incX, double_indexed *Y);
void dmdasum(const int fold, const int N, const double *X, const int incX, double *priY, const int incpriY, double *carY, const int inccarY);
double didssq(const int fold, const int N, const double *X, const int incX, const double scaleY, double_indexed *Y);
double dmdssq(const int fold, const int N, const double *X, const int incX, const double scaleY, double *priY, const int incpriY, double *carY, const int inccarY);
void diddot(const int fold, const int N, const double *X, const int incX, const double *Y, const int incY, double_indexed *Z);
void dmddot(const int fold, const int N, const double *X, const int incX, const double *Y, const int incY, double *manZ, const int incmanZ, double *carZ, const int inccarZ);

void zizsum(const int fold, const int N, const void *X, const int incX, double_indexed *Y);
void zmzsum(const int fold, const int N, const void *X, const int incX, double *priY, const int incpriY, double *carY, const int inccarY);
void dizasum(const int fold, const int N, const void *X, const int incX, double_indexed *Y);
void dmzasum(const int fold, const int N, const void *X, const int incX, double *priY, const int incpriY, double *carY, const int inccarY);
double dizssq(const int fold, const int N, const void *X, const int incX, const double scaleY, double_indexed *Y);
double dmzssq(const int fold, const int N, const void *X, const int incX, const double scaleY, double *priY, const int incpriY, double *carY, const int inccarY);
void zizdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double_indexed *Z);
void zmzdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double *manZ, const int incmanZ, double *carZ, const int inccarZ);
void zizdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double_indexed *Z);
void zmzdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double *manZ, const int incmanZ, double *carZ, const int inccarZ);

void sissum(const int fold, const int N, const float *X, const int incX, float_indexed *Y);
void smssum(const int fold, const int N, const float *X, const int incX, float *priY, const int incpriY, float *carY, const int inccarY);
void sisasum(const int fold, const int N, const float *X, const int incX, float_indexed *Y);
void smsasum(const int fold, const int N, const float *X, const int incX, float *priY, const int incpriY, float *carY, const int inccarY);
float sisssq(const int fold, const int N, const float *X, const int incX, const float scaleY, float_indexed *Y);
float smsssq(const int fold, const int N, const float *X, const int incX, const float scaleY, float *priY, const int incpriY, float *carY, const int inccarY);
void sisdot(const int fold, const int N, const float *X, const int incX, const float *Y, const int incY, float_indexed *Z);
void smsdot(const int fold, const int N, const float *X, const int incX, const float *Y, const int incY, float *manZ, const int incmanZ, float *carZ, const int inccarZ);

void cicsum(const int fold, const int N, const void *X, const int incX, float_indexed *Y);
void cmcsum(const int fold, const int N, const void *X, const int incX, float *priY, const int incpriY, float *carY, const int inccarY);
void sicasum(const int fold, const int N, const void *X, const int incX, float_indexed *Y);
void smcasum(const int fold, const int N, const void *X, const int incX, float *priY, const int incpriY, float *carY, const int inccarY);
float sicssq(const int fold, const int N, const void *X, const int incX, const float scaleY, float_indexed *Y);
float smcssq(const int fold, const int N, const void *X, const int incX, const float scaleY, float *priY, const int incpriY, float *carY, const int inccarY);
void cicdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float_indexed *Z);
void cmcdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float *manZ, const int incmanZ, float *carZ, const int inccarZ);
void cicdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float_indexed *Z);
void cmcdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float *manZ, const int incmanZ, float *carZ, const int inccarZ);

/*

//-==============================-//
//  Accumulator: to improve performance
//  of single aggregation by using
//  Iblas functions
//-==============================-//


typedef struct dIAccumulator {
	Idouble v;
	int counter;
	int size;
	double* BUFFER;
} dIAccum;

typedef struct sIAccumulator {
	Ifloat v;
	int counter;
	int size;
	float* BUFFER;
} sIAccum;

typedef struct zIAccumulator {
	dIcomplex v;
	int counter;
	int size;
	double complex* BUFFER;
} zIAccum;

typedef struct cIAccumulator {
	sIcomplex v;
	int counter;
	int size;
	float complex* BUFFER;
} cIAccum;

extern void dIAccReset(dIAccum *acc);
extern void dIAccInit(dIAccum *acc, int BUFFER_SIZE);
extern void dIAccDestroy(dIAccum *acc);
extern void dIAccumulate(dIAccum *acc, double x);
extern void dIAccumulates(dIAccum *acc, int n, double* x, int inc);
extern double dIAccExtract(dIAccum *acc);

extern void sIAccReset(sIAccum *acc);
extern void sIAccInit(sIAccum *acc, int BUFFER_SIZE);
extern void sIAccDestroy(sIAccum *acc);
extern void sIAccumulate(sIAccum *acc, float x);
extern void sIAccumulates(sIAccum *acc, int n, float* x, int inc);
extern float sIAccExtract(sIAccum *acc);

extern void zIAccReset(zIAccum *acc);
extern void zIAccInit(zIAccum *acc, int BUFFER_SIZE);
extern void zIAccDestroy(zIAccum *acc);
extern void zIAccumulate(zIAccum *acc, double complex x);
extern void zIAccumulates(zIAccum *acc, int n, double complex* x, int inc);
extern double complex zIAccExtract(zIAccum *acc);

extern void cIAccReset(cIAccum *acc);
extern void cIAccInit(cIAccum *acc, int BUFFER_SIZE);
extern void cIAccDestroy(cIAccum *acc);
extern void cIAccumulate(cIAccum *acc, float complex x);
extern void cIAccumulates(cIAccum *acc, int n, float complex* x, int inc);
extern float complex cIAccExtract(cIAccum *acc);
*/

/*
void dgemvI(const int fold, const rblas_order_t Order,
            const rblas_transpose_t TransA, const int M, const int N,
            const double *A, const int lda,
            const double *X, const int incX,
            double_indexed *Y, const int incY);
*/

#endif
