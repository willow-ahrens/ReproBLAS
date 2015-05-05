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

double rdsum(const int N, const double* X, const int incX);
double rdasum(const int N, const double* X, const int incX);
double rdnrm2(const int N, const double* X, const int incX);
double rddot(const int N, const double* X, const int incX, const double* Y, const int incY);

float rsdot(const int N, const float* X, const int incX, const float* Y, const int incY);
float rsasum(const int N, const float* X, const int incX);
float rssum(const int N, const float* X, const int incX);
float rsnrm2(const int N, const float* X, const int incX);

void rzsum_sub(const int N, const void* X, int incX, void *sum);
double rdzasum(const int N, const void* X, const int incX);
double rdznrm2(const int N, const void* X, int incX);
void rzdotc_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotc);
void rzdotu_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotu);

void rcsum_sub(const int N, const void* X, const int incX, void *sum);
float rscasum(const int N, const void* X, const int incX);
float rscnrm2(const int N, const void* X, const int incX);
void rcdotc_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotc);
void rcdotu_sub(const int N, const void* X, const int incX, const void* Y, const int incY, void *dotu);

#endif
