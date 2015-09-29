#include <reproBLAS.h>

#include "../../config.h"

void reproBLAS_sgemv(const char Order, const char TransA,
                     const int M, const int N,
                     const float alpha, const float *A, const int lda,
                     const float *X, const int incX,
                     const float beta, float *Y, const int incY){
  reproBLAS_rsgemv(SIDEFAULTFOLD, Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}
