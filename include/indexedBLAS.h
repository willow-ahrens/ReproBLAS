#ifndef _INDEXED_BLAS__H_
#define _INDEXED_BLAS__H_
#include "indexed.h"
#include "reproBLAS.h"
#include <complex.h>

void dgemvI(int fold, rblas_order_t Order,
            rblas_transpose_t TransA, int M, int N,
            double *A, int lda,
            double *X, int incX,
            I_double *Y, int incY);

float samax(const int N, const float *X, const int incX);
double damax(const int N, const double *X, const int incX);
void camax_sub(const int N, const void *X, const int incX, void *amax);
void zamax_sub(const int N, const void *X, const int incX, void *amax);

float samaxm(const int N, const float *X, const int incX, const float *Y, const int incY);
double damaxm(const int N, const double *X, const int incX, const double *Y, const int incY);
void camaxm_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *amaxm);
void zamaxm_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *amaxm);

void didsum(const int fold, const int N, const double *X, const int incX, double_indexed *Y);
void dmdsum(const int fold, const int N, const double *X, const int incX, double *manY, const int incmanY, double *carY, const int inccarY);
void didasum(const int fold, const int N, const double *X, const int incX, double_indexed *Y);
void dmdasum(const int fold, const int N, const double *X, const int incX, double *manY, const int incmanY, double *carY, const int inccarY);
double didnrm(const int fold, const int N, const double *X, const int incX, double_indexed *Y);
double dmdnrm(const int fold, const int N, const double *X, const int incX, double *manY, const int incmanY, double *carY, const int inccarY);
void diddot(const int fold, const int N, const double *X, const int incX, const double *Y, const int incY, double_indexed *Z);
void dmddot(const int fold, const int N, const double *X, const int incX, const double *Y, const int incY, double *manZ, const int incmanZ, double *carZ, const int inccarZ);

void zizsum(const int fold, const int N, const void *X, const int incX, double_indexed *Y);
void zmzsum(const int fold, const int N, const void *X, const int incX, double *manY, const int incmanY, double *carY, const int inccarY);
void dizasum(const int fold, const int N, const void *X, const int incX, double_indexed *Y);
void dmzasum(const int fold, const int N, const void *X, const int incX, double *manY, const int incmanY, double *carY, const int inccarY);
double diznrm(const int fold, const int N, const void *X, const int incX, double_indexed *Y);
double dmznrm(const int fold, const int N, const void *X, const int incX, double *manY, const int incmanY, double *carY, const int inccarY);
void zizdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double_indexed *Z);
void zmzdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double *manZ, const int incmanZ, double *carZ, const int inccarZ);
void zizdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double_indexed *Z);
void zmzdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double *manZ, const int incmanZ, double *carZ, const int inccarZ);

void sissum(const int fold, const int N, const float *X, const int incX, float_indexed *Y);
void smssum(const int fold, const int N, const float *X, const int incX, float *manY, const int incmanY, float *carY, const int inccarY);
void sisasum(const int fold, const int N, const float *X, const int incX, float_indexed *Y);
void smsasum(const int fold, const int N, const float *X, const int incX, float *manY, const int incmanY, float *carY, const int inccarY);
float sisnrm(const int fold, const int N, const float *X, const int incX, float_indexed *Y);
float smsnrm(const int fold, const int N, const float *X, const int incX, float *manY, const int incmanY, float *carY, const int inccarY);
void sisdot(const int fold, const int N, const float *X, const int incX, const float *Y, const int incY, float_indexed *Z);
void smsdot(const int fold, const int N, const float *X, const int incX, const float *Y, const int incY, float *manZ, const int incmanZ, float *carZ, const int inccarZ);

void cicsum(const int fold, const int N, const void *X, const int incX, float_indexed *Y);
void cmcsum(const int fold, const int N, const void *X, const int incX, float *manY, const int incmanY, float *carY, const int inccarY);
void sicasum(const int fold, const int N, const void *X, const int incX, float_indexed *Y);
void smcasum(const int fold, const int N, const void *X, const int incX, float *manY, const int incmanY, float *carY, const int inccarY);
float sicnrm(const int fold, const int N, const void *X, const int incX, float_indexed *Y);
float smcnrm(const int fold, const int N, const void *X, const int incX, float *manY, const int incmanY, float *carY, const int inccarY);
void cicdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float_indexed *Z);
void cmcdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float *manZ, const int incmanZ, float *carZ, const int inccarZ);
void cicdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float_indexed *Z);
void cmcdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float *manZ, const int incmanZ, float *carZ, const int inccarZ);


#define  dsumI_(N,V,INC,S) dsumI1(N,V,INC,DEFAULT_FOLD,(S).m, (S).c)
#define dasumI_(N,V,INC,S) dasumI1(N,V,INC,DEFAULT_FOLD,(S).m, (S).c)
#define dnrm2I_(N,V,INC,S) dnrm2I1(N,V,INC,DEFAULT_FOLD,(S).m, (S).c)
#define  ddotI_(N,X,INCX,Y,INCY,S) ddotI1(N,X,INCX,Y,INCY,DEFAULT_FOLD,(S).m, (S).c)

extern I_double dsumI (int N, double* v, int inc);
extern I_double dasumI(int N, double* v, int inc);
extern double   dnrm2I(int N, double* v, int inc, I_double* sum);
extern I_double ddotI (int N, double* x, int incx, double* y, int incy);

extern I_float sdotI (int N, float* x, int incx, float* y, int incy);
extern I_float sasumI(int N, float* x, int incx);
extern I_float ssumI (int N, float* x, int incx);
extern float   snrm2I(int N, float* x, int incx, I_float* sum);

extern I_double dzasumI (int N, double complex* v, int inc);
extern I_double_Complex zsumI (int N, double complex* v, int inc);
extern double           dznrm2I (int N, double complex* v, int inc, I_double* sum);
extern I_double_Complex zdotcI(int N, double complex* x, int incx,
                 double complex* y, int incy);
extern I_double_Complex zdotuI(int N, double complex* x, int incx,
                 double complex* y, int incy);

extern I_float         scasumI (int N, float complex* v, int inc);
extern I_float_Complex csumI   (int N, float complex* v, int inc);
extern float           scnrm2I (int N, float complex* v, int inc, I_float* sum);
extern I_float_Complex cdotcI  (int N, float complex* x, int incx, float complex* y, int incy);
extern I_float_Complex cdotuI  (int N, float complex* x, int incx, float complex* y, int incy);

// UNBLOCKED VERSIONS

// ++++ MAX |v[i]| ++++
;
// ---- MAX |v[i]| ----

// ++++ MAX |x[i] * y[i]| ++++
// ---- MAX |x[i] * y[i]| ----


// -- double precision ---
extern void dsumI2 (int n, double* v, int incv,
             int fold, double* sum);
extern void dasumI2(int n, double* v, int incv,
             int fold, double* sum);
extern void dnrm2I2(int n, double* v, int incv, double a,
             int fold, double* sum);
extern void ddotI2 (int n, double* v, int incv, double* y, int incy,
             int fold, double* sum);
// -- single precision ---

extern void ssumI2 (int n, float* v, int incv,
             int fold, float* sum);
extern void sasumI2(int n, float* v, int incv,
             int fold, float* sum);
extern void snrm2I2(int n, float* v, int incv, float a,
             int fold, float* sum);
extern void sdotI2 (int n, float* v, int incv, float* y, int incy,
             int fold, float* sum);
// -- double complex ---

extern void zsumI2 (int n, double complex* v, int incv,
             int fold, double complex* sum);
extern void dzasumI2(int n, double complex* v, int incv,
             int fold, double complex* sum);
extern void dznrm2I2(int n, double complex* v, int incv, double a,
             int fold, double complex* sum);
extern void zdotuI2 (int n, double complex* v, int incv,
             double complex* y, int incy,
             int fold, double complex* sum);
extern void zdotcI2 (int n, double complex* v, int incv,
             double complex* y, int incy,
             int fold, double complex* sum);
// -- single precision ---

extern void csumI2 (int n, float complex* v, int incv,
             int fold, float complex* sum);
extern void scasumI2(int n, float complex* v, int incv,
             int fold, float complex* sum);
extern void scnrm2I2(int n, float complex* v, int incv, float a,
             int fold, float complex* sum);
extern void cdotuI2 (int n, float complex* v, int incv,
             float complex* y, int incy,
             int fold, float complex* sum);
extern void cdotcI2 (int n, float complex* v, int incv,
             float complex* y, int incy,
             int fold, float complex* sum);

// BLOCKING VERSION: LEVEL 1

//++++ DOUBLE +++++
extern void   dsumI1_ (int N, int NB,
              double* v, int inc, int fold, double* sum, double* c);
extern void   dasumI1_(int N, int NB,
              double* v, int inc, int fold, double* sum, double* c);
extern double dnrm2I1_(int N, int NB,
              double* v, int inc, int fold, double* sum, double* c);
extern void   ddotI1_ (int N, int NB, double* x, int incx,
              double* y, int incy, int fold, double* dot, double* c);

#define dsumI1(N,V,INC,K,SUM,C)        dsumI1_(N,1024,V,INC,K,SUM,C)
#define dasumI1(N,V,INC,K,SUM,C)       dasumI1_(N,1024,V,INC,K,SUM,C)
#define dnrm2I1(N,V,INC,K,SUM,C)       dnrm2I1_(N,1024,V,INC,K,SUM,C)
#define ddotI1(N,V,INC,Y,INCY,K,DOT,C) ddotI1_(N,1024,V,INC,Y,INCY,K,DOT,C)
//---- DOUBLE ----

//++++ FLOAT +++++
extern float snrm2I1_(int N, int NB, 
             float* x, int incx, int fold, float* sum, float* c);
extern void  ssumI1_ (int N, int NB,
             float* x, int incx, int fold, float* sum, float* c);
void  sasumI1_(int N, int NB,
             float* x, int incx, int fold, float* sum, float* c);
extern void  sdotI1_ (int N, int NB, float* x, int incx, float* y, int incy,
             int fold, float* dot, float* c);

#define sasumI1(N,V,INC,K,SUM,C)       sasumI1_(N,1024,V,INC,K,SUM,C)
#define ssumI1(N,V,INC,K,SUM,C)        ssumI1_(N,1024,V,INC,K,SUM,C)
#define snrm2I1(N,V,INC,K,SUM,C)       snrm2I1_(N,1024,V,INC,K,SUM,C)
#define sdotI1(N,V,INC,Y,INCY,K,DOT,C) sdotI1_(N,1024,V,INC,Y,INCY,K,DOT,C)
//---- FLOAT -----

//++++ DOUBLE COMPLEX +++++
extern void   dzasumI1_(int N, int NB,
                 double complex* x, int incx, int K, 
				 double* sum, double* c);
extern double dznrm2I1_(int N, int NB,
                 double complex* v, int inc, int fold, 
				 double* sum, double* c);
extern void   zdotI1_  (int N, int NB,
                 double complex* v, int inc, double complex* y, int incy,
                 int fold, 
				 double complex* dot, double complex* c, int conj);
extern void   zsumI1_  (int N, int NB,
                  double complex* v, int inc, int fold, 
			      double complex* sum, double complex* c);

#define dzasumI1(N, X, INCX, K, SUM,C)	\
	if (INCX == 1)	\
		dasumI1(2*N, (double*)X, 1, K, SUM, C);	\
	else	\
		dzasumI1_(N, 1024, X, INCX, K, SUM, C);

#define dznrm2I1(N, X, INCX, K, SUM, C)	\
	((INCX == 1)	? dnrm2I1(2*N, (double*)X, 1, K, SUM, C) : \
		dznrm2I1_(N, 1024, X, INCX, K, SUM, C))

#define zdotcI1(N,X,INCX,Y,INCY,K,DOT,C) \
         zdotI1_(N,1024,X,INCX,Y,INCY,K,DOT,C,1)
#define zdotuI1(N,X,INCX,Y,INCY,K,DOT,C) \
         zdotI1_(N,1024,X,INCX,Y,INCY,K,DOT,C,0)
#define zsumI1(N,X,INCX,K,DOT,C)         zsumI1_(N,1024,X,INCX,K,DOT,C)
//---- DOUBLE COMPLEX ----

//++++ COMPLEX +++++
extern void   scasumI1_(int N, int NB,
                 float complex* x, int inc, int K, 
                 float* sum, float* c, float complex* work);
extern float  scnrm2I1_(int N, int NB,
                 float complex* v, int inc, int K, 
                 float* sum, float* c, float complex* work);
extern void   cdotI1_  (int N, int NB,
                 float complex* v, int inc, float complex* y, int incy,
                 int K,   float complex* dot, float* c, int conj);
extern void   csumI1_  (int N, int NB,
                 float complex* v, int inc, int K, 
                 float complex* sum, float* c);

#define scasumI1(N, X, INCX, K, SUM, C, WORK)	\
	if (INCX == 1)	\
		sasumI1(2*N, (float*)X, 1, K, SUM, C);	\
	else	\
		scasumI1_(N, 1024, X, INCX, K, SUM, C, WORK);

#define scnrm2I1(N, X, INCX, K, SUM, C, WORK)	\
	((INCX == 1)	? snrm2I1_(2*N, 1024, (float*)X, 1, K, SUM, C) : \
		scnrm2I1_(N, 1024, X, INCX, K, SUM, C, WORK))

#define cdotcI1(N,X,INCX,Y,INCY,K,DOT,C) \
         cdotI1_(N,1024,X,INCX,Y,INCY,K,DOT,C,1)
#define cdotuI1(N,X,INCX,Y,INCY,K,DOT,C) \
         cdotI1_(N,1024,X,INCX,Y,INCY,K,DOT,C,0)
#define csumI1(N,X,INCX,K,DOT,C)         csumI1_(N,1024,X,INCX,K,DOT,C)
//---- COMPLEX ----

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

#endif
