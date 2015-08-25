#include <reproBLAS.h>

#include "../../config.h"

void reproBLAS_dgemv(const char Order, const char TransA,
                     const int M, const int N,
                     const double alpha, const double *A, const int lda,
                     const double *X, const int incX,
                     const double beta, double *Y, const int incY){
  reproBLAS_rdgemv(DIDEFAULTFOLD, Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}
