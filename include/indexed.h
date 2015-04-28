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

void cmprint(const int fold, float* repX, int increpX, float* carX, int inccarX);
void ciprint(const int fold, float_complex_indexed *X);
void dmprint(const int fold, double *repX, int increpX, double *carX, int inccarX);
void diprint(const int fold, double_indexed *X);
void smprint(const int fold, float *repX, int increpX, float *carX, int inccarX);
void siprint(const int fold, float_indexed *X);
void zmprint(const int fold, double *repX, int increpX, double *carX, int inccarX);
void ziprint(const int fold, double_complex_indexed *X);

void dmdmadd(const int fold, double *repX, int increpX, double *carX, int inccarX, double* repY, int increpY, double* carY, int inccarY);
void didiadd(const int fold, double_indexed *X, double_indexed *Y);
void zmzmadd(const int fold, double *repX, int increpX, double *carX, int inccarX, double* repY, int increpY, double* carY, int inccarY);
void ziziadd(const int fold, double_complex_indexed *X, double_complex_indexed *Y);
void dmddeposit(const int fold, double X, double *repY, int increpY);
void diddeposit(const int fold, double X, double_indexed *Y);
void dmdadd(const int fold, double X, double *repY, int increpY, double *carY, int inccarY);
void didadd(const int fold, double X, double_indexed *Y);
void zmzdeposit(const int fold, void *X, double *repY, int increpY);
void zizdeposit(const int fold, void *X, double_complex_indexed *Y);
void zmzadd(const int fold, void *X, double *repY, int increpY, double *carY, int inccarY);
void zizadd(const int fold, void *X, double_complex_indexed *Y);

extern int diwidth();
extern int dicapacity();


double dbound(int index);
void dmbound(const int fold, int index, double *repY, int increpY);
int dindex(double X);
int dmindex(double *repX);
int diindex(double_indexed *X);

int smindex(float *repX);
int siindex(float_indexed *X);
int sindex(float X);
double sbound(int index);
void smbound(const int fold, int index, float *repY, int increpY);

void dmdupdate(const int fold, double X, double* repY, int increpY, double* carY, int inccarY);
void didupdate(const int fold, double X, double_indexed *Y);
void zmdupdate(const int fold, double X, double* repY, int increpY, double* carY, int inccarY);
void zidupdate(const int fold, double X, double_complex_indexed *Y);
void zmzupdate(const int fold, void *X, double* repY, int increpY, double* carY, int inccarY);
void zizupdate(const int fold, void *X, double_complex_indexed *Y);
void smsupdate(const int fold, float X, float* repY, int increpY, float* carY, int inccarY);
void sisupdate(const int fold, float X, float_indexed *Y);
void cmsupdate(const int fold, float X, float* repY, int increpY, float* carY, int inccarY);
void cisupdate(const int fold, float X, float_complex_indexed *Y);
void cmcupdate(const int fold, void *X, float* repY, int increpY, float* carY, int inccarY);
void cicupdate(const int fold, void *X, float_complex_indexed *Y);

void dmrenorm(const int fold, double* repX, int increpX, double* carX, int inccarX);
void direnorm(const int fold, double_indexed *X);
void zmrenorm(const int fold, double* repX, int increpX, double* carX, int inccarX);
void zirenorm(const int fold, double_complex_indexed *X);
void smrenorm(const int fold, float* repX, int increpX, float* carX, int inccarX);
void sirenorm(const int fold, float_indexed *X);
void cmrenorm(const int fold, float* repX, int increpX, float* carX, int inccarX);
void cirenorm(const int fold, float_complex_indexed *X);

// COMPUTE THE UNIT IN TH FIRST PLACE
extern double ufp(double x);

// RENORMALIZATION TO AVOID OVERFLOW

// CONVERT A DOUBLE TO INDEXED FORMAT
void didconv(const int fold, double x, double_indexed *y);

// CONVERT AN INDEXED FP BACK TO DOUBLE
double ddiconv(const int fold, double_indexed *x);

// CONVERT AN INDEXED COMPLEX TO COMPLEX
void zziconv_sub(const int fold, double_complex_indexed *x, void *y);

// CONVERT A DOUBLE COMPLEX TO INDEXED FORMAT
extern void zizconv(const int fold, void *x, double_complex_indexed *y);

extern int sicapacity();
extern int siwidth();


void smsmadd(const int fold, float *repX, int increpX, float *carX, int inccarX, float* repY, int increpY, float* carY, int inccarY) ;
void sisiadd(const int fold, float_indexed *X, float_indexed *Y);
void cmcmadd(const int fold, float *repX, int increpX, float *carX, int inccarX, float* repY, int increpY, float* carY, int inccarY) ;
void ciciadd(const int fold, float_complex_indexed *X, float_complex_indexed *Y);
void smsdeposit(const int fold, float X, float *repY, int increpY);
void sisdeposit(const int fold, float X, float_indexed *Y);
void smsadd(const int fold, float X, float *repY, int increpY, float *carY, int inccarY);
void sisadd(const int fold, float X, float_indexed *Y);
void cmcdeposit(const int fold, void *X, float *repY, int increpY);
void cicdeposit(const int fold, void *X, float_complex_indexed *Y);
void cmcadd(const int fold, void *X, float *repY, int increpY, float *carY, int inccarY);
void cicadd(const int fold, void *X, float_complex_indexed *Y);


// NEGATION
extern void sINeg1(int fold, float* x, float* c, int inc);
extern void cINeg1(int fold, float complex* x, float* c, int inc);

// UNIT IN THE FIRST PLACE
extern float ufpf(float x);

// RENORMALIZATION TO AVOID OVERFLOW

// CONVERSION FROM INDEXED FORMAT TO FLOAT
extern float ssiconv(const int fold, float_indexed *x);
extern void cciconv_sub(const int fold, float_complex_indexed *x, void *y);

// CONVERT FROM FLOAT TO INDEXED FORMAT
extern void sisconv(const int fold, float x, float_indexed *y);
extern void cicconv(const int fold, void *x, float_complex_indexed *y);


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
