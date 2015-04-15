#ifndef _REPRODUCIBLE_BLAS__H_
#define _REPRODUCIBLE_BLAS__H_
#define rblas_INDEX size_t  /* this may vary between platforms */
#include <complex.h>

typedef enum rblas_order {
  rblas_Row_Major,
  rblas_Col_Major
} rblas_order_t;
typedef enum rblas_transpose {
  rblas_No_Trans,
  rblas_Trans,
  rblas_Conj_Trans
} rblas_transpose_t;
typedef enum rblas_uplo {
  rblas_Upper,
  rblas_Lower
} rblas_uplo_t;
typedef enum rblas_diag {
  rblas_Non_Unit,
  rblas_Unit
} rblas_diag_t;
typedef enum rblas_side {
  rblas_Left,
  rblas_Right
} rblas_side_t;

void rdgemv(const rblas_order_t order,
            const rblas_transpose_t TransA, const int M, const int N,
            const double *A, const int lda,
            const double *X, const int incX,
            double *Y, const int incY);

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

#endif
