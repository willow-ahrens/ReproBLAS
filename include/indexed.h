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
extern void dIprint1(int n, double *x, double* c, int inc);
extern void zIprint1(int n, double complex* x, double complex* carry, int inc);

extern void dIAdd1(int K,
		double* x, double* xc, int incx,
		double* y, double* yc, int incy);
extern void dIAdd(I_double* X, I_double Y);
extern void zIAdd1(int K,
		double complex* x, double complex* xc, int incx,
		double complex* y, double complex* yc, int incy);
extern void zIAdd(I_double_Complex* X, I_double_Complex Y);

// NEGATION
extern void dINeg1(int fold, double* x, double* c, int inc);
extern void zINeg1(int fold, double complex* x, double complex* c, int inc);

// ADD A DOUBLE TO AN INDEXED FP
// [INPUT]
//    FOLD  : NB OF BINS OF LEADING BITS (NB OF M TO BE COMPUTED)
//    Y     : FP TO BE ADDED TO X
// [IN/OUTPUT]
//    X     : INDEXED FP, AT RETURN X = X + Y
extern void dIAddd1(int fold, double* x, int inc, double y);
extern void dIAddd(I_double* X, double Y);

extern void zIAddz1(int fold, double complex* x, int inc, double complex y);
extern void zIAddz(I_double_Complex* X, double complex Y);

// UPDADATE INDEX FP BY NEW MAXIMUM ABSOLUTE VALUE
extern void dIUpdate1 (int fold, double* x, double* c, int ldx, double y);
extern void zIUpdates1(int fold, double complex* x, double complex* c, int ldx, double y);

extern void zIUpdate1 (int fold, double complex* x, double complex* c, int ldx,double complex y);


double dbound(int index);
void dmbound(int index, double *repY, int increpY, int fold);
int dindex(double X);
int diindex(double_indexed *X);

// COMPUTE THE UNIT IN TH FIRST PLACE
extern double ufp(double x);

// RENORMALIZATION TO AVOID OVERFLOW
// [INPUT]
//    FOLD  : NB OF BINS OF LEADING BITS (NB OF M TO BE COMPUTED)
// [IN/OUTPUT]
//    X     : INDEXED FP
extern void dIRenorm1(int fold, double* X, double* C, int inc);
extern void zIRenorm1(int fold, double complex* rep, double complex* c, int inc);

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

extern void sIprint1(int n, float* x, float* carry, int inc);
extern void cIprint1(int n, float complex* x, float* carry, int inc);

// ADD A FLOAT TO AN INDEXED FP
// [INPUT]
//    FOLD  : NB OF BINS OF LEADING BITS (NB OF M TO BE COMPUTED)
//    Y     : FP TO BE ADDED TO X
// [IN/OUTPUT]
//    X     : INDEXED FP, AT RETURN X = X + Y
extern void sIAddf1(int fold, float* x, int inc, float y);
extern void sIAddf(I_float* X, float Y);
extern void cIAddc1(int fold, float complex* x, int inc, float complex Y);
extern void cIAddc(I_float_Complex* X, float complex Y);

// NEGATION
extern void sINeg1(int fold, float* x, float* c, int inc);
extern void cINeg1(int fold, float complex* x, float* c, int inc);

// COMPUTE BOUNDARIES BASED ON MAXIMUM ABSOLUTE VALUE
// [INPUT]
//    FOLD  : NB OF BINS OF LEADING BITS (NB OF M TO BE COMPUTED)
//    W     : WIDTH OF EACH BINS
//    MAX   : MAXIMUM ABSOLUTE VALUE OF INPUT VALUES TO BE SUMMED
// [OUTPUT]
//    M     : PRECOMPUTED BOUNDARIES
extern int  sIBoundary_(int fold, float max, float* M, int inc);

// UNIT IN THE FIRST PLACE
extern float ufpf(float x);

// X += Y
extern void sIAdd1(int fold, float* x,
	float* xc, int incx, float* y, float* yc, int incy);
extern void sIAdd(I_float* X, I_float Y);

extern void cIAdd1(int fold, float complex* x, float* xc, int incx,
	float complex* y, float* yc, int incy);
extern void cIAdd(I_float_Complex* X, I_float_Complex Y);

//====
extern void sIUpdate1(int fold, float y, float* x, float* C, int inc);

// UPDATE A COMPLEX USING A FLOAT
extern void cIUpdates1(int K, float complex* X, float* C, int INC, float Y);

// UPDATE A COMPLEX USING A COMPLEX
extern void cIUpdate1(int K, float complex* X, float* C,int INC, float complex Y);

// RENORMALIZATION TO AVOID OVERFLOW
// [INPUT]
//    FOLD  : NB OF BINS OF LEADING BITS (NB OF M TO BE COMPUTED)
// [IN/OUTPUT]
//    X     : INDEXED FP
extern void sIRenorm1(int fold, float* X, float* C, int inc);

extern void cIRenorm1(int fold, float complex* sum, float* C, int inc);

// CONVERSION FROM INDEXED FORMAT TO FLOAT
extern float ssiconv(float_indexed *x, int fold);
extern void cciconv_sub(float_complex_indexed *x, void *y, int fold);

// CONVERT FROM FLOAT TO INDEXED FORMAT
extern void sisconv(float x, float_indexed *y, int fold);
extern void cicconv(void *x, float_complex_indexed *y, int fold);


/*******************************************/
/* WRAPPER FOR DOUBLE PRECISION            */
/*******************************************/

#define dIprint(X) dIprint1(DEFAULT_FOLD, (X).m, (X).c, 1)
#define zIprint(X) zIprint1(DEFAULT_FOLD, \
	(double complex*)((X).m), (double complex*)((X).c), 1)

//====================================//
// ADDITION
//====================================//


// ADDING A NATIVE FP TO AN INDEXED FP
#define dIUpdate_(X,Y) dIUpdate1(DEFAULT_FOLD,  (X).m, (X).c, 1, fabs(Y))
#define dIAddd_(X,Y)   dIAddd1(DEFAULT_FOLD, (X).m, 1, Y)
#define dIRenorm_(X)   dIRenorm1(DEFAULT_FOLD, (X).m, (X).c, 1)

#define zIUpdate_(X,Y) zIUpdate1(DEFAULT_FOLD, 	\
		(double complex*)((X).m),(double complex*)((X).c),1,Y)
#define zIAddz_(X,Y) zIAddz1(DEFAULT_FOLD, (double complex*)((X).m), 1, Y)
#define zIRenorm_(X) zIRenorm1(DEFAULT_FOLD,	\
		(double complex*)((X).m), (double complex*)((X).c), 1)

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

#define sIprint(X) sIprint1(DEFAULT_FOLD, (X).m, (X).c, 1)
#define cIprint(X) cIprint1(DEFAULT_FOLD, (float complex*)(X).m, (X).c, 1)

//====================================//
// ADDITION
//====================================//


// ADDING A NATIVE FP TO AN INDEXED FP
#define sIUpdate_(X,Y) sIUpdate1(DEFAULT_FOLD, fabs(Y), (X).m, (X).c, 1)
#define sIAddf_(X,Y)   sIAddf1(DEFAULT_FOLD, (X).m, 1, Y)
#define sIRenorm_(X)   sIRenorm1(DEFAULT_FOLD, (X).m, (X).c, 1)

#define cIUpdate_(X,Y) cIUpdate1(DEFAULT_FOLD,(float complex*)((X).m),(X).c,1,Y)
#define cIAddc_(X,Y) cIAddc1(DEFAULT_FOLD, (float complex*)((X).m), 1, Y)
#define cIRenorm_(X) cIRenorm1(DEFAULT_FOLD,(float complex*)((X).m), (X).c, 1)

//====================================//
// NEGATION
//====================================//
#define sINeg(X) dINeg1(DEFAULT_FOLD, (X).m, (X).c, 1)
#define cINeg(X) zINeg1(DEFAULT_FOLD, (float complex*)(X).m, (X).c, 1)

//====================================//
// CONVERSION
//====================================//

#endif
