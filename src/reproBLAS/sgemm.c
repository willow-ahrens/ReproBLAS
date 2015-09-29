#include <reproBLAS.h>

#include "../../config.h"

void reproBLAS_sgemm(const char Order, const char TransA, const char TransB,
                     const int M, const int N, const int K,
                     const float alpha, const float *A, const int lda,
                     const float *B, const int ldb,
                     const float beta, float *C, const int ldc){
  reproBLAS_rsgemm(SIDEFAULTFOLD, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
