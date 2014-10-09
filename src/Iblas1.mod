#include <complex.h>

// UNBLOCKED VERSIONS

// ++++ MAX |v[i]| ++++
extern float    samax (int n, float* x, int incx);
extern double   damax (int n, double* v, int inc);
extern float complex camax (int n, float complex* x, int incx);
extern double complex zamax (int n, double complex* x, int incx);
// ---- MAX |v[i]| ----

// ++++ MAX |x[i] * y[i]| ++++
extern float  samaxm(int n, float* x, int incx, float* y, int incy);
extern double damaxm(int n, double* x, int incx, double* y, int incy);
extern float  complex camaxm(int n, float complex* x, int incx,
               float complex* y, int incy);
extern double complex zamaxm(int n, double complex* x, int incx,
               double complex* y, int incy);
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
              double* v, int inc, int fold, int W, double* sum, double* c);
extern void   dasumI1_(int N, int NB,
              double* v, int inc, int fold, int W, double* sum, double* c);
extern double dnrm2I1_(int N, int NB,
              double* v, int inc, int fold, int W, double* sum, double* c);
extern void   ddotI1_ (int N, int NB, double* x, int incx,
              double* y, int incy, int fold, int W, double* dot, double* c);

#define dsumI1(N,V,INC,K,W,SUM,C)        dsumI1_(N,1024,V,INC,K,W,SUM,C)
#define dasumI1(N,V,INC,K,W,SUM,C)       dasumI1_(N,1024,V,INC,K,W,SUM,C)
#define dnrm2I1(N,V,INC,K,W,SUM,C)       dnrm2I1_(N,1024,V,INC,K,W,SUM,C)
#define ddotI1(N,V,INC,Y,INCY,K,W,DOT,C) ddotI1_(N,1024,V,INC,Y,INCY,K,W,DOT,C)
//---- DOUBLE ----

//++++ FLOAT +++++
extern float snrm2I1_(int N, int NB, 
             float* x, int incx, int fold, int W, float* sum, F_CARRY_T* c);
extern void  ssumI1_ (int N, int NB,
             float* x, int incx, int fold, int W, float* sum, F_CARRY_T* c);
void  sasumI1_(int N, int NB,
             float* x, int incx, int fold, int W, float* sum, F_CARRY_T* c);
extern void  sdotI1_ (int N, int NB, float* x, int incx, float* y, int incy,
             int fold, int W, float* dot, F_CARRY_T* c);

#define sasumI1(N,V,INC,K,W,SUM,C)       sasumI1_(N,1024,V,INC,K,W,SUM,C)
#define ssumI1(N,V,INC,K,W,SUM,C)        ssumI1_(N,1024,V,INC,K,W,SUM,C)
#define snrm2I1(N,V,INC,K,W,SUM,C)       snrm2I1_(N,1024,V,INC,K,W,SUM,C)
#define sdotI1(N,V,INC,Y,INCY,K,W,DOT,C) sdotI1_(N,1024,V,INC,Y,INCY,K,W,DOT,C)
//---- FLOAT -----

//++++ DOUBLE COMPLEX +++++
extern void   dzasumI1_(int N, int NB,
                 double complex* x, int incx, int K, int W,
				 double* sum, double* c);
extern double dznrm2I1_(int N, int NB,
                 double complex* v, int inc, int fold, int W,
				 double* sum, double* c);
extern void   zdotI1_  (int N, int NB,
                 double complex* v, int inc, double complex* y, int incy,
                 int fold, int W,
				 double complex* dot, double complex* c, int conj);
extern void   zsumI1_  (int N, int NB,
                  double complex* v, int inc, int fold, int W,
			      double complex* sum, double complex* c);

#define dzasumI1(N, X, INCX, K, W, SUM,C)	\
	if (INCX == 1)	\
		dasumI1(2*N, (double*)X, 1, K, W, SUM, C);	\
	else	\
		dzasumI1_(N, 1024, X, INCX, K, W, SUM, C);

#define dznrm2I1(N, X, INCX, K, W, SUM, C)	\
	((INCX == 1)	? dnrm2I1(2*N, (double*)X, 1, K, W, SUM, C) : \
		dznrm2I1_(N, 1024, X, INCX, K, W, SUM, C))

#define zdotcI1(N,X,INCX,Y,INCY,K,W,DOT,C) \
         zdotI1_(N,1024,X,INCX,Y,INCY,K,W,DOT,C,1)
#define zdotuI1(N,X,INCX,Y,INCY,K,W,DOT,C) \
         zdotI1_(N,1024,X,INCX,Y,INCY,K,W,DOT,C,0)
#define zsumI1(N,X,INCX,K,W,DOT,C)         zsumI1_(N,1024,X,INCX,K,W,DOT,C)
//---- DOUBLE COMPLEX ----

//++++ COMPLEX +++++
extern void   scasumI1_(int N, int NB,
                 float complex* x, int inc, int K, int W,
                 float* sum, F_CARRY_T* c, float complex* work);
extern float  scnrm2I1_(int N, int NB,
                 float complex* v, int inc, int K, int W,
                 float* sum, F_CARRY_T* c, float complex* work);
extern void   cdotI1_  (int N, int NB,
                 float complex* v, int inc, float complex* y, int incy,
                 int K, int W,  float complex* dot, F_CARRY_T* c, int conj);
extern void   csumI1_  (int N, int NB,
                 float complex* v, int inc, int K, int W,
                 float complex* sum, F_CARRY_T* c);

#define scasumI1(N, X, INCX, K, W, SUM, C, WORK)	\
	if (INCX == 1)	\
		sasumI1(2*N, (float*)X, 1, K, W, SUM, C);	\
	else	\
		scasumI1_(N, 1024, X, INCX, K, W, SUM, C, WORK);

#define scnrm2I1(N, X, INCX, K, W, SUM, C, WORK)	\
	((INCX == 1)	? snrm2I1_(2*N, 1024, (float*)X, 1, K, W, SUM, C) : \
		scnrm2I1_(N, 1024, X, INCX, K, W, SUM, C, WORK))

#define cdotcI1(N,X,INCX,Y,INCY,K,W,DOT,C) \
         cdotI1_(N,1024,X,INCX,Y,INCY,K,W,DOT,C,1)
#define cdotuI1(N,X,INCX,Y,INCY,K,W,DOT,C) \
         cdotI1_(N,1024,X,INCX,Y,INCY,K,W,DOT,C,0)
#define csumI1(N,X,INCX,K,W,DOT,C)         csumI1_(N,1024,X,INCX,K,W,DOT,C)
//---- COMPLEX ----

