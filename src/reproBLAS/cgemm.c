#include <reproBLAS.h>

#include "../../config.h"

void reproBLAS_cgemm(const char Order, const char TransA, const char TransB,
                     const int M, const int N, const int K,
                     const void *alpha, const void *A, const int lda,
                     const void *B, const int ldb,
                     const void *beta, void *C, const int ldc){
  reproBLAS_rcgemm(SIDEFAULTFOLD, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
