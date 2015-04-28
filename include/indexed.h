//@todo provide a working explanation for indexed types
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
 * Power users of the library may find themselves wanting to manually specify the underlying rep and carry vectors of an indexed type themselves. If you do not know what these are, don't worry about the manually specified indexed types.
 *
 * The goal of using indexed types is to obtain either more accurate or reproducible summation of floating point numbers. Indexed types are composed of several adjacent bins...
 *
 * The parameter @c fold describes how many bins are used in the indexed types supplied to a subroutine. The maximum value for this parameter can be set in config.h. If you are unsure of what value to use for @fold, we recommend 3. Note that the @c fold of indexed types must be the same for all indexed types that interact with each other. Operations on more than one indexed type assume all indexed types being operated upon have the same @c fold. Note that the @c fold of an indexed type may not be changed once the type has been allocated. A common use case would be to set the value of @c fold as a global macro in your code and supply it to all indexed functions that you use.
 *
 * Power users of the library may find themselves wanting to 
 */

#ifndef INDEXED_H
#define INDEXED_H
#include <complex.h>
#include <math.h>
#include <stddef.h>

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

#define DEFAULT_FOLD 3

typedef struct Idouble_ {
	double m[DEFAULT_FOLD];
	double c[DEFAULT_FOLD];
} I_double;

#define Idouble I_double
#define Ifloat  I_float
#define dIcomplex I_double_Complex
#define sIcomplex I_float_Complex

typedef struct dIComplex_ {
	double m[2*DEFAULT_FOLD];
	double c[2*DEFAULT_FOLD];
} I_double_Complex;

typedef struct Ifloat_ {
	float  m[DEFAULT_FOLD];
	float c[DEFAULT_FOLD];
} I_float;

typedef struct sIComplex_{
	float m[2*DEFAULT_FOLD];
	float c[2*DEFAULT_FOLD];
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

int diwidth();
int siwidth();

int dicapacity();
int sicapacity();

double dbound(int index);
int dindex(double X);
int diindex(double_indexed *X);
int dmindex(double *repX);

float sbound(int index);
int sindex(float X);
int siindex(float_indexed *X);
int smindex(float *repX);

void ciprint(const int fold, float_complex_indexed *X);
void cmprint(const int fold, float* repX, int increpX, float* carX, int inccarX);
void diprint(const int fold, double_indexed *X);
void dmprint(const int fold, double *repX, int increpX, double *carX, int inccarX);
void siprint(const int fold, float_indexed *X);
void smprint(const int fold, float *repX, int increpX, float *carX, int inccarX);
void ziprint(const int fold, double_complex_indexed *X);
void zmprint(const int fold, double *repX, int increpX, double *carX, int inccarX);

void didiset(const int fold, const double_indexed *X, double_indexed *Y);
void dmdmset(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double *repY, const int increpY, double *carY, const int inccarY);
void ziziset(const int fold, const double_complex_indexed *X, double_complex_indexed *Y);
void zmzmset(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double *repY, const int increpY, double *carY, const int inccarY);
void zidiset(const int fold, const double_indexed *X, double_complex_indexed *Y);
void zmdmset(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double *repY, const int increpY, double *carY, const int inccarY);
void sisiset(const int fold, const float_indexed *X, float_indexed *Y);
void smsmset(const int fold, const float *repX, const int increpX, const float *carX, const int inccarX, float *repY, const int increpY, float *carY, const int inccarY);
void ciciset(const int fold, const float_complex_indexed *X, float_complex_indexed *Y);
void cmcmset(const int fold, const float *repX, const int increpX, const float *carX, const int inccarX, float *repY, const int increpY, float *carY, const int inccarY);
void cisiset(const int fold, const float_indexed *X, float_complex_indexed *Y);
void cmsmset(const int fold, const float *repX, const int increpX, const float *carX, const int inccarX, float *repY, const int increpY, float *carY, const int inccarY);

void disetzero(const int fold, double_indexed *X);
void dmsetzero(const int fold, double *repX, const int increpX, double *carX, const int inccarX);
void zisetzero(const int fold, double_complex_indexed *X);
void zmsetzero(const int fold, double *repX, const int increpX, double *carX, const int inccarX);
void sisetzero(const int fold, float_indexed *X);
void smsetzero(const int fold, float *repX, const int increpX, float *carX, const int inccarX);
void cisetzero(const int fold, float_complex_indexed *X);
void cmsetzero(const int fold, float *repX, const int increpX, float *carX, const int inccarX);

void didiadd(const int fold, const double_indexed *X, double_indexed *Y);
void dmdmadd(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double* repY, const int increpY, double* carY, const int inccarY);
void ziziadd(const int fold, const double_complex_indexed *X, double_complex_indexed *Y);
void zmzmadd(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double* repY, const int increpY, double* carY, const int inccarY);
void sisiadd(const int fold, float_indexed *X, float_indexed *Y);
void smsmadd(const int fold, float *repX, int increpX, float *carX, int inccarX, float* repY, int increpY, float* carY, int inccarY) ;
void ciciadd(const int fold, float_complex_indexed *X, float_complex_indexed *Y);
void cmcmadd(const int fold, float *repX, int increpX, float *carX, int inccarX, float* repY, int increpY, float* carY, int inccarY) ;

void didadd(const int fold, const double X, double_indexed *Y);
void dmdadd(const int fold, const double X, double *repY, const int increpY, double *carY, const int inccarY);
void zizadd(const int fold, const void *X, double_complex_indexed *Y);
void zmzadd(const int fold, const void *X, double *repY, const int increpY, double *carY, const int inccarY);
void sisadd(const int fold, float X, float_indexed *Y);
void smsadd(const int fold, float X, float *repY, int increpY, float *carY, int inccarY);
void cicadd(const int fold, void *X, float_complex_indexed *Y);
void cmcadd(const int fold, void *X, float *repY, int increpY, float *carY, int inccarY);

void didupdate(const int fold, double X, double_indexed *Y);
void dmdupdate(const int fold, double X, double* repY, int increpY, double* carY, int inccarY);
void zizupdate(const int fold, void *X, double_complex_indexed *Y);
void zmzupdate(const int fold, void *X, double* repY, int increpY, double* carY, int inccarY);
void zidupdate(const int fold, double X, double_complex_indexed *Y);
void zmdupdate(const int fold, double X, double* repY, int increpY, double* carY, int inccarY);
void sisupdate(const int fold, float X, float_indexed *Y);
void smsupdate(const int fold, float X, float* repY, int increpY, float* carY, int inccarY);
void cicupdate(const int fold, void *X, float_complex_indexed *Y);
void cmcupdate(const int fold, void *X, float* repY, int increpY, float* carY, int inccarY);
void cisupdate(const int fold, float X, float_complex_indexed *Y);
void cmsupdate(const int fold, float X, float* repY, int increpY, float* carY, int inccarY);

void diddeposit(const int fold, const double X, double_indexed *Y);
void dmddeposit(const int fold, const double X, double *repY, const int increpY);
void zizdeposit(const int fold, const void *X, double_complex_indexed *Y);
void zmzdeposit(const int fold, const void *X, double *repY, const int increpY);
void sisdeposit(const int fold, float X, float_indexed *Y);
void smsdeposit(const int fold, float X, float *repY, int increpY);
void cicdeposit(const int fold, void *X, float_complex_indexed *Y);
void cmcdeposit(const int fold, void *X, float *repY, int increpY);

void dmrenorm(const int fold, double* repX, int increpX, double* carX, int inccarX);
void direnorm(const int fold, double_indexed *X);
void zmrenorm(const int fold, double* repX, int increpX, double* carX, int inccarX);
void zirenorm(const int fold, double_complex_indexed *X);
void smrenorm(const int fold, float* repX, int increpX, float* carX, int inccarX);
void sirenorm(const int fold, float_indexed *X);
void cmrenorm(const int fold, float* repX, int increpX, float* carX, int inccarX);
void cirenorm(const int fold, float_complex_indexed *X);


void didconv(const int fold, double x, double_indexed *y);
void dmdconv(const int fold, double x, double* repy, int increpy, double* cary, int inccary);
void zizconv(const int fold, void *x, double_complex_indexed *y);
void zmzconv(const int fold, void *x, double *repy, int increpy, double *cary, int inccary);
void sisconv(const int fold, float x, float_indexed *y);
void smsconv(const int fold, float x, float* repy, int increpy, float* cary, int inccary);
void cicconv(const int fold, void *x, float_complex_indexed *y);
void cmcconv(const int fold, void *x, float *repy, int increpy, float *cary, int inccary);
double ddiconv(const int fold, double_indexed *x);
double ddmconv(const int fold, double* repx, int increpx, double* carx, int inccarx);
void zziconv_sub(const int fold, double_complex_indexed *x, void *y);
void zzmconv_sub(const int fold, double *repx, int increpx, double *carx, int inccarx, void *y);
float ssiconv(const int fold, float_indexed *x);
float ssmconv(const int fold, float* repx, int increpx, float* carx, int inccarx);
void ccmconv_sub(const int fold, float *repx, int increpx, float *carx, int inccarx, void *y);
void cciconv_sub(const int fold, float_complex_indexed *x, void *y);

void dmnegate(const int fold, double* repX, const int increpX, double* carX, const int inccarX);
void dinegate(const int fold, double_indexed* X);
void zmnegate(const int fold, double* repX, const int increpX, double* carX, const int inccarX);
void zinegate(const int fold, double_complex_indexed* X);
void smnegate(const int fold, float* repX, const int increpX, float* carX, const int inccarX);
void sinegate(const int fold, float_indexed* X);
void cmnegate(const int fold, float* repX, const int increpX, float* carX, const int inccarX);
void cinegate(const int fold, float_complex_indexed* X);

double ufp(double x);
float ufpf(float x);
void dmbound(const int fold, int index, double *repY, int increpY);
void smbound(const int fold, int index, float *repY, int increpY);
#endif
