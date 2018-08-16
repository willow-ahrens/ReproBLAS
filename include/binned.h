/**
 * @file  binned.h
 * @brief binned.h defines the binned types and the lower level functions associated with their use.
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
 * The goal of using binned types is to obtain either more accurate or reproducible summation of floating point numbers. In reproducible summation, floating point numbers are split into several slices along predefined boundaries in the exponent range. The space between two boundaries is called a bin. Binned types are composed of several accumulators, each accumulating the slices in a particular bin. The accumulators correspond to the largest consecutive nonzero bins seen so far.
 *
 * The parameter @c fold describes how many accumulators are used in the binned types supplied to a subroutine (an binned type with @c k accumulators  is @c k-fold). The default value for this parameter can be set in config.h. If you are unsure of what value to use for @c fold, we recommend 3. Note that the @c fold of binned types must be the same for all binned types that interact with each other. Operations on more than one binned type assume all binned types being operated upon have the same @c fold. Note that the @c fold of an binned type may not be changed once the type has been allocated. A common use case would be to set the value of @c fold as a global macro in your code and supply it to all binned functions that you use.
 *
 * @internal
 * Power users of the library may find themselves wanting to manually specify the underlying primary and carry vectors of an binned type themselves. If you do not know what these are, don't worry about the manually specified binned types.
 */

#ifndef BINNED_H_
#define BINNED_H_
#include <stddef.h>
#include <stdlib.h>
#include <float.h>

/**
 * @brief The binned double datatype
 *
 * To allocate a #double_binned, call binned_dballoc()
 *
 * @warning A #double_binned is, under the hood, an array of @c double. Therefore, if you have defined an array of #double_binned, you must index it by multiplying the index into the array by the number of underlying @c double that make up the #double_binned. This number can be obtained by a call to binned_dbnum()
 */
typedef double double_binned;

/**
 * @brief The binned complex double datatype
 *
 * To allocate a #double_complex_binned, call binned_zballoc()
 *
 * @warning A #double_complex_binned is, under the hood, an array of @c double. Therefore, if you have defined an array of #double_complex_binned, you must index it by multiplying the index into the array by the number of underlying @c double that make up the #double_complex_binned. This number can be obtained by a call to binned_zbnum()
 */
typedef double double_complex_binned;

/**
 * @brief The binned float datatype
 *
 * To allocate a #float_binned, call binned_sballoc()
 *
 * @warning A #float_binned is, under the hood, an array of @c float. Therefore, if you have defined an array of #float_binned, you must index it by multiplying the index into the array by the number of underlying @c float that make up the #float_binned. This number can be obtained by a call to binned_sbnum()
 */
typedef float float_binned;

/**
 * @brief The binned complex float datatype
 *
 * To allocate a #float_complex_binned, call binned_cballoc()
 *
 * @warning A #float_complex_binned is, under the hood, an array of @c float. Therefore, if you have defined an array of #float_complex_binned, you must index it by multiplying the index into the array by the number of underlying @c float that make up the #float_complex_binned. This number can be obtained by a call to binned_cbnum()
 */
typedef float float_complex_binned;

/**
 * @brief Binned double precision bin width
 *
 * bin width (in bits)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
#define DBWIDTH 40

/**
 * @brief Binned single precision bin width
 *
 * bin width (in bits)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
#define SBWIDTH 13

/**
 * @brief Binned double precision maximum index
 *
 * maximum index (inclusive)
 *
 * @author Peter Ahrens
 * @date   24 Jun 2015
 */
#define binned_DBMAXINDEX (((DBL_MAX_EXP - DBL_MIN_EXP + DBL_MANT_DIG - 1)/DBWIDTH) - 1)

/**
 * @brief Binned single precision maximum index
 *
 * maximum index (inclusive)
 *
 * @author Peter Ahrens
 * @date   24 Jun 2015
 */
#define binned_SBMAXINDEX (((FLT_MAX_EXP - FLT_MIN_EXP + FLT_MANT_DIG - 1)/SBWIDTH) - 1)

/**
 * @brief The maximum double precision fold supported by the library.
 *
 * @author Peter Ahrens
 * @date   14 Jan 2016
 */
#define binned_DBMAXFOLD (binned_DBMAXINDEX + 1)

/**
 * @brief The maximum single precision fold supported by the library.
 *
 * @author Peter Ahrens
 * @date   14 Jan 2016
 */
#define binned_SBMAXFOLD (binned_SBMAXINDEX + 1)

/**
 * @brief Binned double precision deposit endurance
 *
 * The number of deposits that can be performed before a renorm is necessary. Applies also to binned complex double precision.
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
#define binned_DBENDURANCE (1 << (DBL_MANT_DIG - DBWIDTH - 2))

/**
 * @brief Binned single precision deposit endurance
 *
 * The number of deposits that can be performed before a renorm is necessary. Applies also to binned complex single precision.
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
#define binned_SBENDURANCE (1 << (FLT_MANT_DIG - SBWIDTH - 2))

/**
 * @brief Binned double precision capacity
 *
 * The maximum number of double precision numbers that can be summed using binned double precision. Applies also to binned complex double precision.
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
#define binned_DBCAPACITY (binned_DBENDURANCE*(1.0/DBL_EPSILON - 1.0))

/**
 * @brief Binned single precision capacity
 *
 * The maximum number of single precision numbers that can be summed using binned single precision. Applies also to binned complex double precision.
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
#define binned_SBCAPACITY (binned_SBENDURANCE*(1.0/FLT_EPSILON - 1.0))

/**
 * @internal
 * @brief Binned double precision compression factor
 *
 * This factor is used to scale down inputs before deposition into the bin of highest index
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
#define binned_DMCOMPRESSION (1.0/(1 << (DBL_MANT_DIG - DBWIDTH + 1)))

/**
 * @internal
 * @brief Binned single precision compression factor
 *
 * This factor is used to scale down inputs before deposition into the bin of highest index
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
#define binned_SMCOMPRESSION (1.0/(1 << (FLT_MANT_DIG - SBWIDTH + 1)))

/**
 * @internal
 * @brief Binned double precision expansion factor
 *
 * This factor is used to scale up inputs after deposition into the bin of highest index
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
#define binned_DMEXPANSION (1.0*(1 << (DBL_MANT_DIG - DBWIDTH + 1)))

/**
 * @internal
 * @brief Binned single precision expansion factor
 *
 * This factor is used to scale up inputs after deposition into the bin of highest index
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
#define binned_SMEXPANSION (1.0*(1 << (FLT_MANT_DIG - SBWIDTH + 1)))

size_t binned_dbsize(const int fold);
size_t binned_zbsize(const int fold);
size_t binned_sbsbze(const int fold);
size_t binned_cbsize(const int fold);

double_binned *binned_dballoc(const int fold);
double_complex_binned *binned_zballoc(const int fold);
float_binned *binned_sballoc(const int fold);
float_complex_binned *binned_cballoc(const int fold);

int binned_dbnum(const int fold);
int binned_zbnum(const int fold);
int binned_sbnum(const int fold);
int binned_cbnum(const int fold);

double binned_dbbound(const int fold, const int N, const double X, const double S);
float binned_sbbound(const int fold, const int N, const float X, const float S);

const double *binned_dmbins(const int X);
const float *binned_smbins(const int X);

int binned_dindex(const double X);
int binned_dmindex(const double *priX);
int binned_dmindex0(const double *priX);

int binned_sindex(const float X);
int binned_smindex(const float *priX);
int binned_smindex0(const float *priX);

int binned_dmdenorm(const int fold, const double *priX);
int binned_zmdenorm(const int fold, const double *priX);
int binned_smdenorm(const int fold, const float *priX);
int binned_cmdenorm(const int fold, const float *priX);

void binned_dbprint(const int fold, const double_binned *X);
void binned_dmprint(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX);
void binned_zbprint(const int fold, const double_complex_binned *X);
void binned_zmprint(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX);
void binned_sbprint(const int fold, const float_binned *X);
void binned_smprint(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX);
void binned_cbprint(const int fold, const float_complex_binned *X);
void binned_cmprint(const int fold, const float* priX, const int incpriX, const float* carX, const int inccarX);

void binned_dbdbset(const int fold, const double_binned *X, double_binned *Y);
void binned_dmdmset(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, double *priY, const int incpriY, double *carY, const int inccarY);
void binned_zbzbset(const int fold, const double_complex_binned *X, double_complex_binned *Y);
void binned_zmzmset(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, double *priY, const int incpriY, double *carY, const int inccarY);
void binned_zbdbset(const int fold, const double_binned *X, double_complex_binned *Y);
void binned_zmdmset(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, double *priY, const int incpriY, double *carY, const int inccarY);
void binned_sbsbset(const int fold, const float_binned *X, float_binned *Y);
void binned_smsmset(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, float *priY, const int incpriY, float *carY, const int inccarY);
void binned_cbcbset(const int fold, const float_complex_binned *X, float_complex_binned *Y);
void binned_cmcmset(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, float *priY, const int incpriY, float *carY, const int inccarY);
void binned_cbsbset(const int fold, const float_binned *X, float_complex_binned *Y);
void binned_cmsmset(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, float *priY, const int incpriY, float *carY, const int inccarY);

void binned_dbsetzero(const int fold, double_binned *X);
void binned_dmsetzero(const int fold, double *priX, const int incpriX, double *carX, const int inccarX);
void binned_zbsetzero(const int fold, double_complex_binned *X);
void binned_zmsetzero(const int fold, double *priX, const int incpriX, double *carX, const int inccarX);
void binned_sbsetzero(const int fold, float_binned *X);
void binned_smsetzero(const int fold, float *priX, const int incpriX, float *carX, const int inccarX);
void binned_cbsetzero(const int fold, float_complex_binned *X);
void binned_cmsetzero(const int fold, float *priX, const int incpriX, float *carX, const int inccarX);

void binned_dbdbadd(const int fold, const double_binned *X, double_binned *Y);
void binned_dmdmadd(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, double* priY, const int incpriY, double* carY, const int inccarY);
void binned_zbzbadd(const int fold, const double_complex_binned *X, double_complex_binned *Y);
void binned_zmzmadd(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, double* priY, const int incpriY, double* carY, const int inccarY);
void binned_sbsbadd(const int fold, const float_binned *X, float_binned *Y);
void binned_smsmadd(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, float* priY, const int incpriY, float* carY, const int inccarY) ;
void binned_cbcbadd(const int fold, const float_complex_binned *X, float_complex_binned *Y);
void binned_cmcmadd(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, float* priY, const int incpriY, float* carY, const int inccarY) ;

void binned_dbdbaddv(const int fold, const int N, const double_binned *X, const int incX, double_binned *Y, const int incY);
void binned_zbzbaddv(const int fold, const int N, const double_complex_binned *X, const int incX, double_complex_binned *Y, const int incY);
void binned_sbsbaddv(const int fold, const int N, const float_binned *X, const int incX, float_binned *Y, const int incY);
void binned_cbcbaddv(const int fold, const int N, const float_complex_binned *X, const int incX, float_complex_binned *Y, const int incY);

void binned_dbdadd(const int fold, const double X, double_binned *Y);
void binned_dmdadd(const int fold, const double X, double *priY, const int incpriY, double *carY, const int inccarY);
void binned_zbzadd(const int fold, const void *X, double_complex_binned *Y);
void binned_zmzadd(const int fold, const void *X, double *priY, const int incpriY, double *carY, const int inccarY);
void binned_sbsadd(const int fold, const float X, float_binned *Y);
void binned_smsadd(const int fold, const float X, float *priY, const int incpriY, float *carY, const int inccarY);
void binned_cbcadd(const int fold, const void *X, float_complex_binned *Y);
void binned_cmcadd(const int fold, const void *X, float *priY, const int incpriY, float *carY, const int inccarY);

void binned_dbdupdate(const int fold, const double X, double_binned *Y);
void binned_dmdupdate(const int fold, const double X, double* priY, const int incpriY, double* carY, const int inccarY);
void binned_zbzupdate(const int fold, const void *X, double_complex_binned *Y);
void binned_zmzupdate(const int fold, const void *X, double* priY, const int incpriY, double* carY, const int inccarY);
void binned_zbdupdate(const int fold, const double X, double_complex_binned *Y);
void binned_zmdupdate(const int fold, const double X, double* priY, const int incpriY, double* carY, const int inccarY);
void binned_sbsupdate(const int fold, const float X, float_binned *Y);
void binned_smsupdate(const int fold, const float X, float* priY, const int incpriY, float* carY, const int inccarY);
void binned_cbcupdate(const int fold, const void *X, float_complex_binned *Y);
void binned_cmcupdate(const int fold, const void *X, float* priY, const int incpriY, float* carY, const int inccarY);
void binned_cbsupdate(const int fold, const float X, float_complex_binned *Y);
void binned_cmsupdate(const int fold, const float X, float* priY, const int incpriY, float* carY, const int inccarY);

void binned_dbddeposit(const int fold, const double X, double_binned *Y);
void binned_dmddeposit(const int fold, const double X, double *priY, const int incpriY);
void binned_zbzdeposit(const int fold, const void *X, double_complex_binned *Y);
void binned_zmzdeposit(const int fold, const void *X, double *priY, const int incpriY);
void binned_sbsdeposit(const int fold, const float X, float_binned *Y);
void binned_smsdeposit(const int fold, const float X, float *priY, const int incpriY);
void binned_cbcdeposit(const int fold, const void *X, float_complex_binned *Y);
void binned_cmcdeposit(const int fold, const void *X, float *priY, const int incpriY);

void binned_dbrenorm(const int fold, double_binned *X);
void binned_dmrenorm(const int fold, double* priX, const int incpriX, double* carX, const int inccarX);
void binned_zbrenorm(const int fold, double_complex_binned *X);
void binned_zmrenorm(const int fold, double* priX, const int incpriX, double* carX, const int inccarX);
void binned_sbrenorm(const int fold, float_binned *X);
void binned_smrenorm(const int fold, float* priX, const int incpriX, float* carX, const int inccarX);
void binned_cbrenorm(const int fold, float_complex_binned *X);
void binned_cmrenorm(const int fold, float* priX, const int incpriX, float* carX, const int inccarX);

void binned_dbdconv(const int fold, const double X, double_binned *Y);
void binned_dmdconv(const int fold, const double X, double* priY, const int incpriY, double* carY, const int inccarY);
void binned_zbzconv(const int fold, const void *X, double_complex_binned *Y);
void binned_zmzconv(const int fold, const void *X, double *priY, const int incpriY, double *carY, const int inccarY);
void binned_sbsconv(const int fold, const float X, float_binned *Y);
void binned_smsconv(const int fold, const float X, float* priY, const int incpriY, float* carY, const int inccarY);
void binned_cbcconv(const int fold, const void *X, float_complex_binned *Y);
void binned_cmcconv(const int fold, const void *X, float *priY, const int incpriY, float *carY, const int inccarY);
double binned_ddbconv(const int fold, const double_binned *X);
double binned_ddmconv(const int fold, const double* priX, const int incpriX, const double* carX, const int inccarX);
void binned_zzbconv_sub(const int fold, const double_complex_binned *X, void *conv);
void binned_zzmconv_sub(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, void *conv);
float binned_ssbconv(const int fold, const float_binned *X);
float binned_ssmconv(const int fold, const float* priX, const int incpriX, const float* carX, const int inccarX);
void binned_ccbconv_sub(const int fold, const float_complex_binned *X, void *conv);
void binned_ccmconv_sub(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, void *conv);

void binned_dbnegate(const int fold, double_binned* X);
void binned_dmnegate(const int fold, double* priX, const int incpriX, double* carX, const int inccarX);
void binned_zbnegate(const int fold, double_complex_binned* X);
void binned_zmnegate(const int fold, double* priX, const int incpriX, double* carX, const int inccarX);
void binned_sbnegate(const int fold, float_binned* X);
void binned_smnegate(const int fold, float* priX, const int incpriX, float* carX, const int inccarX);
void binned_cbnegate(const int fold, float_complex_binned* X);
void binned_cmnegate(const int fold, float* priX, const int incpriX, float* carX, const int inccarX);

double binned_dscale(const double X);
float binned_sscale(const float X);

void binned_dmdrescale(const int fold, const double X, const double scaleY, double *priY, const int incpriY, double *carY, const int inccarY);
void binned_zmdrescale(const int fold, const double X, const double scaleY, double *priY, const int incpriY, double *carY, const int inccarY);
void binned_smsrescale(const int fold, const float X, const float scaleY, float *priY, const int incpriY, float *carY, const int inccarY);
void binned_cmsrescale(const int fold, const float X, const float scaleY, float *priY, const int incpriY, float *carY, const int inccarY);

double binned_dmdmaddsq(const int fold, const double scaleX, const double *priX, const int incpriX, const double *carX, const int inccarX, const double scaleY, double *priY, const int incpriY, double *carY, const int inccarY);
double binned_dbdbaddsq(const int fold, const double scaleX, const double_binned *X, const double scaleY, double_binned *Y);
float binned_smsmaddsq(const int fold, const float scaleX, const float *priX, const int incpriX, const float *carX, const int inccarX, const float scaleY, float *priY, const int incpriY, float *carY, const int inccarY);
float binned_sbsbaddsq(const int fold, const float scaleX, const float_binned *X, const float scaleY, float_binned *Y);

double binned_ufp(const double X);
float binned_ufpf(const float X);
#endif
