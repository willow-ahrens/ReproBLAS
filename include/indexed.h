#ifndef INDEXED_H
#define INDEXED_H
#include <complex.h>
#include <math.h>
#include <stddef.h>

typedef double double_indexed;
typedef double double_complex_indexed;
typedef float float_indexed;
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
int dinum(const int fold);
int zinum(const int fold);
double_indexed *dialloc(const int fold);
double_complex_indexed *zialloc(const int fold);
void dmdmset(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double *repY, const int increpY, double *carY, const int inccarY);
void didiset(const int fold, const double_indexed *X, double_indexed *Y);
void zmzmset(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double *repY, const int increpY, double *carY, const int inccarY);
void ziziset(const int fold, const double_complex_indexed *X, double_complex_indexed *Y);
void dmzmset(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double *repY, const int increpY, double *carY, const int inccarY);
void diziset(const int fold, const double_indexed *X, double_complex_indexed *Y);
void dmsetzero(const int fold, double *repX, const int increpX, double *carX, const int inccarX);
void disetzero(const int fold, double_indexed *X);
void zmsetzero(const int fold, double *repX, const int increpX, double *carX, const int inccarX);
void zisetzero(const int fold, double_complex_indexed *X);
size_t sisize(const int fold);
size_t cisize(const int fold);
int sinum(const int fold);
int cinum(const int fold);
float_indexed *sialloc(const int fold);
float_complex_indexed *cialloc(const int fold);
void smsmset(const int fold, const float *repX, const int increpX, const float *carX, const int inccarX, float *repY, const int increpY, float *carY, const int inccarY);
void sisiset(const int fold, const float_indexed *X, float_indexed *Y);
void cmcmset(const int fold, const float *repX, const int increpX, const float *carX, const int inccarX, float *repY, const int increpY, float *carY, const int inccarY);
void ciciset(const int fold, const float_complex_indexed *X, float_complex_indexed *Y);
void cmsmset(const int fold, const float *repX, const int increpX, const float *carX, const int inccarX, float *repY, const int increpY, float *carY, const int inccarY);
void cisiset(const int fold, const float_indexed *X, float_complex_indexed *Y);
void smsetzero(const int fold, float *repX, const int increpX, float *carX, const int inccarX);
void sisetzero(const int fold, float_indexed *X);
void cmsetzero(const int fold, float *repX, const int increpX, float *carX, const int inccarX);
void cisetzero(const int fold, float_complex_indexed *X);

// PRINT

void cmprint(float* repX, int increpX, float* carX, int inccarX, int fold);
void ciprint(float_complex_indexed *X, int fold);
void dmprint(double *repX, int increpX, double *carX, int inccarX, int fold);
void diprint(double_indexed *X, int fold);
void smprint(float *repX, int increpX, float *carX, int inccarX, int fold);
void siprint(float_indexed *X, int fold);
void zmprint(double *repX, int increpX, double *carX, int inccarX, int fold);
void ziprint(double_complex_indexed *X, int fold);

void dmdmadd(double *repX, int increpX, double *carX, int inccarX, double* repY, int increpY, double* carY, int inccarY, int fold);
void didiadd(double_indexed *X, double_indexed *Y, int fold);
void zmzmadd(double *repX, int increpX, double *carX, int inccarX, double* repY, int increpY, double* carY, int inccarY, int fold);
void ziziadd(double_complex_indexed *X, double_complex_indexed *Y, int fold);
void dmddeposit(double X, double *repY, int increpY, int fold);
void diddeposit(double X, double_indexed *Y, int fold);
void dmdadd(double X, double *repY, int increpY, double *carY, int inccarY, int fold);
void didadd(double X, double_indexed *Y, int fold);
void zmzdeposit(void *X, double *repY, int increpY, int fold);
void zizdeposit(void *X, double_complex_indexed *Y, int fold);
void zmzadd(void *X, double *repY, int increpY, double *carY, int inccarY, int fold);
void zizadd(void *X, double_complex_indexed *Y, int fold);

extern int diwidth();
extern int dicapacity();


double dbound(int index);
void dmbound(int index, double *repY, int increpY, int fold);
int dindex(double X);
int dmindex(double *repX);
int diindex(double_indexed *X);

int smindex(float *repX);
int siindex(float_indexed *X);
int sindex(float X);
double sbound(int index);
void smbound(int index, float *repY, int increpY, int fold);

void dmdupdate(double X, double* repY, int increpY, double* carY, int inccarY, int fold);
void didupdate(double X, double_indexed *Y, int fold);
void zmdupdate(double X, double* repY, int increpY, double* carY, int inccarY, int fold);
void zidupdate(double X, double_complex_indexed *Y, int fold);
void zmzupdate(void *X, double* repY, int increpY, double* carY, int inccarY, int fold);
void zizupdate(void *X, double_complex_indexed *Y, int fold);
void smsupdate(float X, float* repY, int increpY, float* carY, int inccarY, int fold);
void sisupdate(float X, float_indexed *Y, int fold);
void cmsupdate(float X, float* repY, int increpY, float* carY, int inccarY, int fold);
void cisupdate(float X, float_complex_indexed *Y, int fold);
void cmcupdate(void *X, float* repY, int increpY, float* carY, int inccarY, int fold);
void cicupdate(void *X, float_complex_indexed *Y, int fold);

void dmrenorm(double* repX, int increpX, double* carX, int inccarX, int fold);
void direnorm(double_indexed *X, int fold);
void zmrenorm(double* repX, int increpX, double* carX, int inccarX, int fold);
void zirenorm(double_complex_indexed *X, int fold);
void smrenorm(float* repX, int increpX, float* carX, int inccarX, int fold);
void sirenorm(float_indexed *X, int fold);
void cmrenorm(float* repX, int increpX, float* carX, int inccarX, int fold);
void cirenorm(float_complex_indexed *X, int fold);

// COMPUTE THE UNIT IN TH FIRST PLACE
extern double ufp(double x);

// RENORMALIZATION TO AVOID OVERFLOW

// CONVERT A DOUBLE TO INDEXED FORMAT
void didconv(double x, double_indexed *y, int fold);

// CONVERT AN INDEXED FP BACK TO DOUBLE
double ddiconv(double_indexed *x, int fold);

// CONVERT AN INDEXED COMPLEX TO COMPLEX
void zziconv_sub(double_complex_indexed *x, void *y, int fold);

// CONVERT A DOUBLE COMPLEX TO INDEXED FORMAT
extern void zizconv(void *x, double_complex_indexed *y, int fold);

extern int sicapacity();
extern int siwidth();


void smsmadd(float *repX, int increpX, float *carX, int inccarX, float* repY, int increpY, float* carY, int inccarY, int fold) ;
void sisiadd(float_indexed *X, float_indexed *Y, int fold);
void cmcmadd(float *repX, int increpX, float *carX, int inccarX, float* repY, int increpY, float* carY, int inccarY, int fold) ;
void ciciadd(float_complex_indexed *X, float_complex_indexed *Y, int fold);
void smsdeposit(float X, float *repY, int increpY, int fold);
void sisdeposit(float X, float_indexed *Y, int fold);
void smsadd(float X, float *repY, int increpY, float *carY, int inccarY, int fold);
void sisadd(float X, float_indexed *Y, int fold);
void cmcdeposit(void *X, float *repY, int increpY, int fold);
void cicdeposit(void *X, float_complex_indexed *Y, int fold);
void cmcadd(void *X, float *repY, int increpY, float *carY, int inccarY, int fold);
void cicadd(void *X, float_complex_indexed *Y, int fold);


// NEGATION
extern void sINeg1(int fold, float* x, float* c, int inc);
extern void cINeg1(int fold, float complex* x, float* c, int inc);

// UNIT IN THE FIRST PLACE
extern float ufpf(float x);

// RENORMALIZATION TO AVOID OVERFLOW

// CONVERSION FROM INDEXED FORMAT TO FLOAT
extern float ssiconv(float_indexed *x, int fold);
extern void cciconv_sub(float_complex_indexed *x, void *y, int fold);

// CONVERT FROM FLOAT TO INDEXED FORMAT
extern void sisconv(float x, float_indexed *y, int fold);
extern void cicconv(void *x, float_complex_indexed *y, int fold);


/*******************************************/
/* WRAPPER FOR DOUBLE PRECISION            */
/*******************************************/


//====================================//
// ADDITION
//====================================//



//====================================//
// NEGATION
//====================================//

//====================================//
// CONVERSION
//====================================//


/*******************************************/
/* WRAPPER FOR SINGLE PRECISION            */
/*******************************************/


//====================================//
// ADDITION
//====================================//




//====================================//
// NEGATION
//====================================//

//====================================//
// CONVERSION
//====================================//


void dmnegate(const int fold, double* repX, const int increpX, double* carX, const int inccarX);
void dinegate(const int fold, double_indexed* X);
void zmnegate(const int fold, double* repX, const int increpX, double* carX, const int inccarX);
void zinegate(const int fold, double_complex_indexed* X);
void smnegate(const int fold, float* repX, const int increpX, float* carX, const int inccarX);
void sinegate(const int fold, float_indexed* X);
void cmnegate(const int fold, float* repX, const int increpX, float* carX, const int inccarX);
void cinegate(const int fold, float_complex_indexed* X);
#endif
