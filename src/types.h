#ifndef _REPRODUCIBLE_BLAS_TYPES__H_
#define _REPRODUCIBLE_BLAS_TYPES__H_

#define _MATH__COMPLEX_

// DOUBLE COMPLEX
#ifdef _MATH__COMPLEX_
#include <complex.h>

typedef double complex dcomplex;
#define ZREAL_(X) creal(X)
#define ZIMAG_(X) cimag(X)
#define ZSET_(X, Real_, Imag_) X = Real_ + Imag_ * _Complex_I;
#define ZSET_REAL_(X, Real_) X = Real_ + ZIMAG_(X) * _Complex_I;
#define ZSET_IMAG_(X, Imag_) X = ZREAL_(X) + Imag_ * _Complex_I;
#define ZSCAL_(Y,X,a) Y = a * X;

#undef I

#else
///////////////////////////////////////
// DATA TYPE
///////////////////////////////////////

typedef struct dcomplex_ {
	double real;
	double imag;
} dcomplex;

#define ZREAL_(X) X.real
#define ZIMAG_(X) X.imag
#define ZSET_(X, R, I) { X.real = R; X.imag = I; }
#define ZSET_REAL_(X, R) X.real = R;
#define ZSET_IMAG_(X, R) X.imag = I;
#define ZSCAL_(Y,X,a) { Y.real = a * X.real; Y.imag = a * X.imag; }

#endif

// DOUBLE COMPLEX
#ifdef _MATH__COMPLEX_
#include <complex.h>
#undef I

typedef float  complex scomplex;
#define CREAL_(X) crealf(X)
#define CIMAG_(X) cimagf(X)
#define CSET_(X, Real_, Imag_) X = Real_ + Imag_ * _Complex_I;
#define CSET_REAL_(X, Real_) X = Real_ + CIMAG_(X) * _Complex_I;
#define CSET_IMAG_(X, Imag_) X = CREAL_(X) + Imag_ * _Complex_I;
#define CSCAL_(Y,X,a) Y = a * X;

#else
///////////////////////////////////////
// DATA TYPE
///////////////////////////////////////

typedef struct scomplex_ {
	float real;
	float imag;
} scomplex;
#define CREAL_(X) X.real
#define CIMAG_(X) X.imag
#define CSET_(X, R, I) { X.real = R; X.imag = I; }
#define CSET_REAL_(X, R) X.real = R;
#define CSET_IMAG_(X, I) X.imag = I;
#define CSCAL_(Y,X,a) { Y.real = a * X.real; Y.imag = a * X.imag; }

#endif

#endif
