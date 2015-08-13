#include <reproBLAS.h>

#include "../../config.h"

void reproBLAS_zgemv(const char Order,
                     const char TransA, const int M, const int N,
                     const void *alpha, const void *A, const int lda,
                     const void *X, const int incX,
                     const void *beta, void *Y, const int incY){
  reproBLAS_rzgemv(DIDEFAULTFOLD, Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}
