#include <reproBLAS.h>

#include "../../config.h"

void reproBLAS_dgemm(const char Order, const char TransA, const char TransB,
                     const int M, const int N, const int K,
                     const double alpha, const double *A, const int lda,
                     const double *B, const int ldb,
                     const double beta, double *C, const int ldc){
  reproBLAS_rdgemm(DIDEFAULTFOLD, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
