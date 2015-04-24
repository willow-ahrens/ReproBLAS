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

static inline size_t dISize(int fold){
  return 2*fold*sizeof(double);
}

static inline size_t zISize(int fold){
  return 4*fold*sizeof(double);
}

static inline size_t sISize(int fold){
  return 2*fold*sizeof(float);
}

static inline size_t cISize(int fold){
  return 4*fold*sizeof(float);
}

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

// SET ZERO
#define ISetZero_(K,M,C) {	\
	int i;	\
	for (i = 0; i < K; i++) { \
		M[i] = 0;	\
		C[i] = 0;	\
	}	\
}
#define dISetZero(I) ISetZero_(DEFAULT_FOLD, (I).m, (I).c)
#define sISetZero(I) ISetZero_(DEFAULT_FOLD, (I).m, (I).c)
#define zISetZero(I) ISetZero_(2*DEFAULT_FOLD, (I).m, (I).c)
#define cISetZero(I) ISetZero_(2*DEFAULT_FOLD, (I).m, (I).c)

// DI = SI
#define ISet_(K,DI,SI) {	\
	int i;	\
	for (i = 0; i <  K; i++) {	\
		(DI).m[i] = (SI).m[i]; \
		(DI).c[i] = (SI).c[i]; \
	}	\
}
#define dISet(DI,SI) ISet_(DEFAULT_FOLD, DI, SI)
#define sISet(DI,SI) ISet_(DEFAULT_FOLD, DI, SI)
#define zISet(DI,SI) ISet_(2*DEFAULT_FOLD, DI, SI)
#define cISet(DI,SI) ISet_(2*DEFAULT_FOLD, DI, SI)

// ZI = DI
#define zdISet_(K,ZI,DI) {	\
	ISetZero_(2 * K, (ZI).m, (ZI).c)\
	int i;	\
	for (i = 0; i <  K; i++) {	\
		(ZI).m[2 * i] = (DI).m[i]; \
		(ZI).c[2 * i] = (DI).c[i]; \
	}	\
}
#define zdISet(ZI, DI) zdISet_(DEFAULT_FOLD, ZI, DI)
#define csISet(CI, SI) zdISet_(DEFAULT_FOLD, CI, SI)

extern int dIWidth();
extern int dICapacity();


// NEGATION
extern void dINeg1(int fold, double* x, double* c, int inc);
extern void zINeg1(int fold, double complex* x, double complex* c, int inc);


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

extern int sICapacity();
extern int sIWidth();


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
#define dINeg(X) dINeg1(DEFAULT_FOLD, (X).m, (X).c, 1)
#define zINeg(X) zINeg1(DEFAULT_FOLD, (double complex*)(X).m, (double complex*)(X).c, 1)

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
#define sINeg(X) dINeg1(DEFAULT_FOLD, (X).m, (X).c, 1)
#define cINeg(X) zINeg1(DEFAULT_FOLD, (float complex*)(X).m, (X).c, 1)

//====================================//
// CONVERSION
//====================================//

#endif
