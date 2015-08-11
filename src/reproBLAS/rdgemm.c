#include <stdlib.h>

#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

void rdgemm(const char Order, const char TransA, const char TransB,
            const int M, const int N, const int K,
            const double alpha, const double *A, const int lda,
            const double *B, const int ldb,
            const double beta, double *C, const int ldc){
  double_indexed *CI;
  int i;
  int j;

  CI = (double_indexed*)malloc(N * M * disize(DIDEFAULTFOLD));
  switch(Order){
    case 'r':
    case 'R':
      for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
          didconv(DIDEFAULTFOLD, C[i * lda + j] * beta, CI + (i * lda + j) * dinum(DIDEFAULTFOLD));
        }
      }
      didgemm(DIDEFAULTFOLD, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, N);
      for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
          C[i * lda + j] = ddiconv(DIDEFAULTFOLD, CI + (i * lda + j) * dinum(DIDEFAULTFOLD));
        }
      }
      break;
    default:
      for(j = 0; j < N; j++){
        for(i = 0; i < M; i++){
          didconv(DIDEFAULTFOLD, C[j * lda + i] * beta, CI + (j * lda + i) * dinum(DIDEFAULTFOLD));
        }
      }
      didgemm(DIDEFAULTFOLD, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, M);
      for(j = 0; j < N; j++){
        for(i = 0; i < M; i++){
          C[j * lda + i] = ddiconv(DIDEFAULTFOLD, CI + (j * lda + i) * dinum(DIDEFAULTFOLD));
        }
      }
      break;
  }
  free(CI);
}
