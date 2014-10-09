#define  dsumI_(N,V,INC,S) dsumI1(N,V,INC,DEFAULT_FOLD,0,(S).m, (S).c)
#define dasumI_(N,V,INC,S) dasumI1(N,V,INC,DEFAULT_FOLD,0,(S).m, (S).c)
#define dnrm2I_(N,V,INC,S) dnrm2I1(N,V,INC,DEFAULT_FOLD,0,(S).m, (S).c)
#define  ddotI_(N,X,INCX,Y,INCY,S) ddotI1(N,X,INCX,Y,INCY,DEFAULT_FOLD,0,(S).m, (S).c)

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

// SEQUENTIAL REPRODUCIBLE VERSIONS
extern double rdsum (int N, double* v, int inc);
extern double rdasum(int N, double* v, int inc);
extern double rdnrm2(int N, double* v, int inc);
extern double rddot (int N, double* x, int incx, double* y, int incy);

extern float  rsdot (int N, float* x, int incx, float* y, int incy);
extern float  rsasum(int N, float* x, int incx);
extern float  rssum  (int N, float* x, int incx);
extern float  rsnrm2(int N, float* x, int incx);

extern double         rdzasum (int N, double complex* v, int inc);
extern double complex rzsum   (int N, double complex* v, int inc);
extern double         rdznrm2 (int N, double complex* v, int inc);
extern double complex rzdotc  (int N, double complex* x, int incx, double complex* y, int incy);
extern double complex rzdotu  (int N, double complex* x, int incx, double complex* y, int incy);

extern float         rscasum (int N, float complex* v, int inc);
extern float complex rcsum   (int N, float complex* v, int inc);
extern float         rscnrm2 (int N, float complex* v, int inc);
extern float complex rcdotc  (int N, float complex* x, int incx, float complex* y, int incy);
extern float complex rcdotu  (int N, float complex* x, int incx, float complex* y, int incy);

