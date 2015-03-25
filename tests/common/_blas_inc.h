#ifndef __BLAS_INC_H
#define __BLAS_INC_H

#include <complex.h>

#ifdef CBLAS
	extern int cblas_idamax(int, double*, int);
	extern int cblas_isamax(int, float*, int);
	extern int cblas_izamax(int, double complex*, int);
	extern int cblas_icamax(int, float complex*, int);

	extern double cblas_dasum  (int, double*, int);
	extern float  cblas_sasum  (int, float*, int);
	extern double cblas_dzasum (int, double complex*, int);
	extern float  cblas_scasum (int, float complex*, int);

	extern double cblas_dnrm2  (int, double*, int);
	extern float  cblas_snrm2  (int, float*, int);
	extern double cblas_dznrm2 (int, double complex*, int);
	extern float  cblas_scnrm2 (int, float complex*, int);

	extern double cblas_ddot (int, double*, int, double*, int);
	extern float cblas_sdot  (int, float*, int, float*, int);
	extern void cblas_zdotc_sub(int, double complex*, int, double complex*, int, double complex*);
	extern void cblas_zdotu_sub(int, double complex*, int, double complex*, int, double complex*);
	extern void cblas_cdotc_sub(int, float complex*, int, float complex*, int, float complex*);
	extern void cblas_cdotu_sub(int, float complex*, int, float complex*, int, float complex*);

	extern void cblas_dgemv (char,int,int,double,double*,int,double*,int,double,double*,int);

#	define CALL_IDAMAX(R, N, V, INC) R = cblas_idamax(N, V, INC)
#	define CALL_ISAMAX(R, N, V, INC) R = cblas_isamax(N, V, INC)
#	define CALL_IZAMAX(R, N, V, INC) R = cblas_izamax(N, V, INC)
#	define CALL_ICAMAX(R, N, V, INC) R = cblas_icamax(N, V, INC)

#	define CALL_DASUM(R, N, V, INC)  R = cblas_dasum (N, V, INC)
#	define CALL_SASUM(R, N, V, INC)  R = cblas_sasum (N, V, INC)
#	define CALL_DZASUM(R, N, V, INC) R = cblas_dzasum (N, V, INC)
#	define CALL_SCASUM(R, N, V, INC) R = cblas_scasum (N, V, INC)

#	define CALL_DNRM2(R, N, V, INC)   R = cblas_dnrm2 (N, V, INC)
#	define CALL_SNRM2(R, N, V, INC)   R = cblas_snrm2 (N, V, INC)
#	define CALL_DZNRM2(R, N, V, INC)  R = cblas_dznrm2 (N, V, INC)
#	define CALL_SCNRM2(R, N, V, INC)  R = cblas_scnrm2 (N, V, INC)

#	define CALL_DDOT(R,N, V, INC, Y, INCY) R = cblas_ddot (N, V, INC, Y, INCY)
#	define CALL_SDOT(R,N, V, INC, Y, INCY) R = cblas_sdot (N, V, INC, Y, INCY)
#	define CALL_ZDOTC(R,N,V,INC,Y,INCY)  cblas_zdotc_sub(N,V,INC,Y,INCY,&R)
#	define CALL_ZDOTU(R,N,V,INC,Y,INCY)  cblas_zdotu_sub(N,V,INC,Y,INCY,&R)
#	define CALL_CDOTC(R,N,V,INC,Y,INCY)  cblas_cdotc_sub(N,V,INC,Y,INCY,&R)
#	define CALL_CDOTU(R,N,V,INC,Y,INCY)  cblas_cdotu_sub(N,V,INC,Y,INCY,&R)

#	define CALL_DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) \
		cblas_dgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY);
#	define CALL_DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) \
		cblas_dgemm(TRANSA,TRANSB,M,N,KALPHA,A,LDA,B,LDB,BETA,C,LDC);

#elif defined (BLAS)

	extern int idamax_(int*, double*, int*);
	extern int isamax_(int*, float*, int*);
	extern int izamax_(int*, double complex*, int*);
	extern int icamax_(int*, float complex*, int*);

	extern double dasum_  (int*, double*, int*);
	extern float  sasum_  (int*, float*, int*);
	extern double dzasum_ (int*, double complex*, int*);
	extern float  scasum_ (int*, float complex*, int*);

	extern double dnrm2_  (int*, double*, int*);
	extern float  snrm2_  (int*, float* , int*);
	extern double dznrm2_ (int*, double complex*, int*);
	extern float  scnrm2_ (int*, float complex*, int*);

	extern double ddot_ (int*, double*, int*, double*, int*);
	extern float  sdot_ (int*, float*, int*, float*, int*);
	extern double complex zdotc_ (int*, double complex*, int*,
		double complex*, int*);
	extern double complex zdotu_ (int*, double complex*, int*,
		double complex*, int*);
	extern float  complex cdotc_ (int*, float  complex*, int*,
		float  complex*, int*);
	extern float  complex cdotu_ (int*, float  complex*, int*,
		float  complex*, int*);

	extern double dgemv_ (char*,int*,int*,double*,double*,int*,
		double*,int*,double*,double*,int*);
	extern double dgemm_ (char*,char*,int*,int*,int*,double*,
		double*,int*,double*,int*,double*,double*,int*);


#	define CALL_IDAMAX(R, N, V, INC) R = idamax_(&N, V, &INC)
#	define CALL_ISAMAX(R, N, V, INC) R = isamax_(&N, V, &INC)
#	define CALL_IZAMAX(R, N, V, INC) R = izamax_(&N, V, &INC)
#	define CALL_ICAMAX(R, N, V, INC) R = icamax_(&N, V, &INC)

#	define CALL_DASUM(R, N, V, INC)  R = dasum_ (&N, V, &INC)
#	define CALL_DZASUM(R, N, V, INC) R = dzasum_ (&N, V, &INC)
#	define CALL_SASUM(R, N, V, INC)  R = sasum_ (&N, V, &INC)
#	define CALL_SCASUM(R, N, V, INC) R = scasum_ (&N, V, &INC)

#	define CALL_DNRM2(R, N, V, INC)  R = dnrm2_ (&N, V, &INC)
#	define CALL_DZNRM2(R, N, V, INC) R = dznrm2_ (&N, V, &INC)
#	define CALL_SNRM2(R, N, V, INC)  R = snrm2_ (&N, V, &INC)
#	define CALL_SCNRM2(R, N, V, INC) R = scnrm2_ (&N, V, &INC)

#	define CALL_DDOT(R, N, V, INC, Y, INCY)  R = ddot_ (&N, V, &INC, Y, &INCY)
#	define CALL_ZDOTC(R, N, V, INC, Y, INCY) R = zdotc_ (&N, V, &INC, Y, &INCY)
#	define CALL_ZDOTU(R, N, V, INC, Y, INCY) R = zdotu_ (&N, V, &INC, Y, &INCY)
#	define CALL_SDOT(R, N, V, INC, Y, INCY)  R = sdot_ (&N, V, &INC, Y, &INCY)
#	define CALL_CDOTC(R, N, V, INC, Y, INCY) R = cdotc_ (&N, V, &INC, Y, &INCY)
#	define CALL_CDOTU(R, N, V, INC, Y, INCY) R = cdotu_ (&N, V, &INC, Y, &INCY)

#	define CALL_DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) \
		dgemv_(&TRANS,&M,&N,&ALPHA,A,&LDA,X,&INCX,&BETA,Y,&INCY);
#	define CALL_DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) \
		dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);

#endif 

#ifdef SCALAPACK
    extern void blacs_pinfo_( int*, int* );
    extern void blacs_barrier_(int*, char*);
    extern void blacs_get_( int*, int*, int* );
    extern void blacs_gridinit_( int*, char*, int*, int* );
    extern void pdgemm_( char*, char*, int*, int*, int*, double*, double*, int*, int*, int*, double*, int*, int*, int*, double*, double*, int*, int*, int* );  
#endif
#endif
