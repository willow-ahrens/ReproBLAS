//TODO provide a working explanation for indexed types
/**
 * @file  indexed.h
 * @brief indexed.h defines the indexed types and the lower level functions associated with their use.
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

#ifndef INDEXED_H
#define INDEXED_H
#include <stddef.h>
#include <stdlib.h>
#include <float.h>

/**
 * @brief The indexed double datatype
 *
 * To allocate a #double_indexed, call dialloc()
 *
 * @warning A #double_indexed is, under the hood, an array of @c double. Therefore, if you have defined an array of #double_indexed, you must index it by multiplying the index into the array by the number of underlying @c double that make up the #double_indexed. This number can be obtained by a call to dinum()
 */
typedef double double_indexed;

/**
 * @brief The indexed complex double datatype
 *
 * To allocate a #double_complex_indexed, call zialloc()
 *
 * @warning A #double_complex_indexed is, under the hood, an array of @c double. Therefore, if you have defined an array of #double_complex_indexed, you must index it by multiplying the index into the array by the number of underlying @c double that make up the #double_complex_indexed. This number can be obtained by a call to zinum()
 */
typedef double double_complex_indexed;

/**
 * @brief The indexed float datatype
 *
 * To allocate a #float_indexed, call sialloc()
 *
 * @warning A #float_indexed is, under the hood, an array of @c float. Therefore, if you have defined an array of #float_indexed, you must index it by multiplying the index into the array by the number of underlying @c float that make up the #float_indexed. This number can be obtained by a call to sinum()
 */
typedef float float_indexed;

/**
 * @brief The indexed complex float datatype
 *
 * To allocate a #float_complex_indexed, call cialloc()
 *
 * @warning A #float_complex_indexed is, under the hood, an array of @c float. Therefore, if you have defined an array of #float_complex_indexed, you must index it by multiplying the index into the array by the number of underlying @c float that make up the #float_complex_indexed. This number can be obtained by a call to cinum()
 */
typedef float float_complex_indexed;

/**
 * @brief Indexed double precision bin width
 *
 * bin width (in bits)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
#define DIWIDTH 41

/**
 * @brief Indexed single precision bin width
 *
 * bin width (in bits)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
#define SIWIDTH 13

/**
 * @brief Indexed double precision maximum index
 *
 * maximum index (inclusive)
 *
 * @author Peter Ahrens
 * @date   24 Jun 2015
 */
#define DIMAXINDEX (((((DBL_MAX_EXP - DBL_MIN_EXP) + DBL_MANT_DIG) - 2)/DIWIDTH) - 1)

/**
 * @brief Indexed single precision maximum index
 *
 * maximum index (inclusive)
 *
 * @author Peter Ahrens
 * @date   24 Jun 2015
 */
#define SIMAXINDEX (((((FLT_MAX_EXP - FLT_MIN_EXP) + FLT_MANT_DIG) - 2)/SIWIDTH) - 1)

/**
 * @brief Indexed double precision deposit endurance
 *
 * The number of deposits that can be performed before a renorm is necessary. Applies also to indexed complex double precision.
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
#define DIENDURANCE (1 << (DBL_MANT_DIG - DIWIDTH - 2))

/**
 * @brief Indexed single precision deposit endurance
 *
 * The number of deposits that can be performed before a renorm is necessary. Applies also to indexed complex single precision.
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
#define SIENDURANCE (1 << (FLT_MANT_DIG - SIWIDTH - 2))

/**
 * @internal
 * @brief Indexed double precision compression factor
 *
 * This factor is used to scale down inputs before deposition into the bin of highest index
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
#define DMCOMPRESSION (1.0/(1 << (DBL_MANT_DIG - DIWIDTH + 1)))

/**
 * @internal
 * @brief Indexed single precision compression factor
 *
 * This factor is used to scale down inputs before deposition into the bin of highest index
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
#define SMCOMPRESSION (1.0/(1 << (FLT_MANT_DIG - SIWIDTH + 1)))

/**
 * @internal
 * @brief Indexed double precision expansion factor
 *
 * This factor is used to scale up inputs after deposition into the bin of highest index
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
#define DMEXPANSION (1.0*(1 << (DBL_MANT_DIG - DIWIDTH + 1)))

/**
 * @internal
 * @brief Indexed single precision expansion factor
 *
 * This factor is used to scale up inputs after deposition into the bin of highest index
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
#define SMEXPANSION (1.0*(1 << (FLT_MANT_DIG - SIWIDTH + 1)))

/**
 * @internal
 */
typedef struct Idouble_ {
	double m[3];
	double c[3];
} I_double;

/**
 * @internal
 */
#define Idouble I_double
/**
 * @internal
 */
#define Ifloat  I_float
/**
 * @internal
 */
#define dIcomplex I_double_Complex
/**
 * @internal
 */
#define sIcomplex I_float_Complex

/**
 * @internal
 */
typedef struct dIComplex_ {
	double m[2*3];
	double c[2*3];
} I_double_Complex;

/**
 * @internal
 */
typedef struct Ifloat_ {
	float  m[3];
	float c[3];
} I_float;

/**
 * @internal
 */
typedef struct sIComplex_{
	float m[2*3];
	float c[2*3];
} I_float_Complex;

size_t disize(const int fold);
size_t zisize(const int fold);
size_t sisize(const int fold);
size_t cisize(const int fold);

double_indexed *dialloc(const int fold);
double_complex_indexed *zialloc(const int fold);
float_indexed *sialloc(const int fold);
float_complex_indexed *cialloc(const int fold);

int dinum(const int fold);
int zinum(const int fold);
int sinum(const int fold);
int cinum(const int fold);

double dibound(const int fold, const int N, const double X);
float sibound(const int fold, const int N, const float X);

const double *dmbins(const int X);
const float *smbins(const int X);

int dindex(const double X);
int dmindex(const double *priX);
int dmindex0(const double *priX);

int sindex(const float X);
int smindex(const float *priX);
int smindex0(const float *priX);

int dmdenorm(const int fold, const double *priX);
int zmdenorm(const int fold, const double *priX);
int smdenorm(const int fold, const float *priX);
int cmdenorm(const int fold, const float *priX);

void diprint(const int fold, const double_indexed *X);
void dmprint(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX);
void ziprint(const int fold, const double_complex_indexed *X);
void zmprint(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX);
void siprint(const int fold, const float_indexed *X);
void smprint(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX);
void ciprint(const int fold, const float_complex_indexed *X);
void cmprint(const int fold, const float* priX, const int incpriX, const float* carX, const int inccarX);

void didiset(const int fold, const double_indexed *X, double_indexed *Y);
void dmdmset(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, double *priY, const int incpriY, double *carY, const int inccarY);
void ziziset(const int fold, const double_complex_indexed *X, double_complex_indexed *Y);
void zmzmset(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, double *priY, const int incpriY, double *carY, const int inccarY);
void zidiset(const int fold, const double_indexed *X, double_complex_indexed *Y);
void zmdmset(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, double *priY, const int incpriY, double *carY, const int inccarY);
void sisiset(const int fold, const float_indexed *X, float_indexed *Y);
void smsmset(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, float *priY, const int incpriY, float *carY, const int inccarY);
void ciciset(const int fold, const float_complex_indexed *X, float_complex_indexed *Y);
void cmcmset(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, float *priY, const int incpriY, float *carY, const int inccarY);
void cisiset(const int fold, const float_indexed *X, float_complex_indexed *Y);
void cmsmset(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, float *priY, const int incpriY, float *carY, const int inccarY);

void disetzero(const int fold, double_indexed *X);
void dmsetzero(const int fold, double *priX, const int incpriX, double *carX, const int inccarX);
void zisetzero(const int fold, double_complex_indexed *X);
void zmsetzero(const int fold, double *priX, const int incpriX, double *carX, const int inccarX);
void sisetzero(const int fold, float_indexed *X);
void smsetzero(const int fold, float *priX, const int incpriX, float *carX, const int inccarX);
void cisetzero(const int fold, float_complex_indexed *X);
void cmsetzero(const int fold, float *priX, const int incpriX, float *carX, const int inccarX);

void didiadd(const int fold, const double_indexed *X, double_indexed *Y);
void dmdmadd(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, double* priY, const int incpriY, double* carY, const int inccarY);
void ziziadd(const int fold, const double_complex_indexed *X, double_complex_indexed *Y);
void zmzmadd(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, double* priY, const int incpriY, double* carY, const int inccarY);
void sisiadd(const int fold, const float_indexed *X, float_indexed *Y);
void smsmadd(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, float* priY, const int incpriY, float* carY, const int inccarY) ;
void ciciadd(const int fold, const float_complex_indexed *X, float_complex_indexed *Y);
void cmcmadd(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, float* priY, const int incpriY, float* carY, const int inccarY) ;

void didadd(const int fold, const double X, double_indexed *Y);
void dmdadd(const int fold, const double X, double *priY, const int incpriY, double *carY, const int inccarY);
void zizadd(const int fold, const void *X, double_complex_indexed *Y);
void zmzadd(const int fold, const void *X, double *priY, const int incpriY, double *carY, const int inccarY);
void sisadd(const int fold, const float X, float_indexed *Y);
void smsadd(const int fold, const float X, float *priY, const int incpriY, float *carY, const int inccarY);
void cicadd(const int fold, const void *X, float_complex_indexed *Y);
void cmcadd(const int fold, const void *X, float *priY, const int incpriY, float *carY, const int inccarY);

void didupdate(const int fold, const double X, double_indexed *Y);
void dmdupdate(const int fold, const double X, double* priY, const int incpriY, double* carY, const int inccarY);
void zizupdate(const int fold, const void *X, double_complex_indexed *Y);
void zmzupdate(const int fold, const void *X, double* priY, const int incpriY, double* carY, const int inccarY);
void zidupdate(const int fold, const double X, double_complex_indexed *Y);
void zmdupdate(const int fold, const double X, double* priY, const int incpriY, double* carY, const int inccarY);
void sisupdate(const int fold, const float X, float_indexed *Y);
void smsupdate(const int fold, const float X, float* priY, const int incpriY, float* carY, const int inccarY);
void cicupdate(const int fold, const void *X, float_complex_indexed *Y);
void cmcupdate(const int fold, const void *X, float* priY, const int incpriY, float* carY, const int inccarY);
void cisupdate(const int fold, const float X, float_complex_indexed *Y);
void cmsupdate(const int fold, const float X, float* priY, const int incpriY, float* carY, const int inccarY);

void diddeposit(const int fold, const double X, double_indexed *Y);
void dmddeposit(const int fold, const double X, double *priY, const int incpriY);
void zizdeposit(const int fold, const void *X, double_complex_indexed *Y);
void zmzdeposit(const int fold, const void *X, double *priY, const int incpriY);
void sisdeposit(const int fold, const float X, float_indexed *Y);
void smsdeposit(const int fold, const float X, float *priY, const int incpriY);
void cicdeposit(const int fold, const void *X, float_complex_indexed *Y);
void cmcdeposit(const int fold, const void *X, float *priY, const int incpriY);

void direnorm(const int fold, double_indexed *X);
void dmrenorm(const int fold, double* priX, const int incpriX, double* carX, const int inccarX);
void zirenorm(const int fold, double_complex_indexed *X);
void zmrenorm(const int fold, double* priX, const int incpriX, double* carX, const int inccarX);
void sirenorm(const int fold, float_indexed *X);
void smrenorm(const int fold, float* priX, const int incpriX, float* carX, const int inccarX);
void cirenorm(const int fold, float_complex_indexed *X);
void cmrenorm(const int fold, float* priX, const int incpriX, float* carX, const int inccarX);

void didconv(const int fold, const double X, double_indexed *Y);
void dmdconv(const int fold, const double X, double* priY, const int incpriY, double* carY, const int inccarY);
void zizconv(const int fold, const void *X, double_complex_indexed *Y);
void zmzconv(const int fold, const void *X, double *priY, const int incpriY, double *carY, const int inccarY);
void sisconv(const int fold, const float X, float_indexed *Y);
void smsconv(const int fold, const float X, float* priY, const int incpriY, float* carY, const int inccarY);
void cicconv(const int fold, const void *X, float_complex_indexed *Y);
void cmcconv(const int fold, const void *X, float *priY, const int incpriY, float *carY, const int inccarY);
double ddiconv(const int fold, const double_indexed *X);
double ddmconv(const int fold, const double* priX, const int incpriX, const double* carX, const int inccarX);
void zziconv_sub(const int fold, const double_complex_indexed *X, void *conv);
void zzmconv_sub(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, void *conv);
float ssiconv(const int fold, const float_indexed *X);
float ssmconv(const int fold, const float* priX, const int incpriX, const float* carX, const int inccarX);
void cciconv_sub(const int fold, const float_complex_indexed *X, void *conv);
void ccmconv_sub(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, void *conv);

void dinegate(const int fold, double_indexed* X);
void dmnegate(const int fold, double* priX, const int incpriX, double* carX, const int inccarX);
void zinegate(const int fold, double_complex_indexed* X);
void zmnegate(const int fold, double* priX, const int incpriX, double* carX, const int inccarX);
void sinegate(const int fold, float_indexed* X);
void smnegate(const int fold, float* priX, const int incpriX, float* carX, const int inccarX);
void cinegate(const int fold, float_complex_indexed* X);
void cmnegate(const int fold, float* priX, const int incpriX, float* carX, const int inccarX);

double dscale(const double X);
float sscale(const float X);

void dmdrescale(const int fold, const double X, const double scaleY, double *priY, const int incpriY, double *carY, const int inccarY);
void zmdrescale(const int fold, const double X, const double scaleY, double *priY, const int incpriY, double *carY, const int inccarY);
void smsrescale(const int fold, const float X, const float scaleY, float *priY, const int incpriY, float *carY, const int inccarY);
void cmsrescale(const int fold, const float X, const float scaleY, float *priY, const int incpriY, float *carY, const int inccarY);

double dmdmaddsq(const int fold, const double scaleX, const double *priX, const int incpriX, const double *carX, const int inccarX, const double scaleY, double* priY, const int incpriY, double* carY, const int inccarY);
float smsmaddsq(const int fold, const float scaleX, const float *priX, const int incpriX, const float *carX, const int inccarX, const float scaleY, float* priY, const int incpriY, float* carY, const int inccarY);

double ufp(const double X);
float ufpf(const float X);
#endif
