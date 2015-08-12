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

  if(M == 0 || N == 0){
    return;
  }

  CI = (double_indexed*)malloc(M * N * disize(DIDEFAULTFOLD));
  switch(Order){
    case 'r':
    case 'R':
      if(beta == 1.0){
        for(i = 0; i < M; i++){
          for(j = 0; j < N; j++){
            didconv(DIDEFAULTFOLD, C[i * ldc + j], CI + (i * N + j) * dinum(DIDEFAULTFOLD));
          }
        }
      }else if(beta == 0.0){
        memset(CI, 0.0, M * N * disize(DIDEFAULTFOLD));
      }else{
        for(i = 0; i < M; i++){
          for(j = 0; j < N; j++){
            didconv(DIDEFAULTFOLD, C[i * ldc + j] * beta, CI + (i * N + j) * dinum(DIDEFAULTFOLD));
          }
        }
      }
      didgemm(DIDEFAULTFOLD, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, N);
      for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
          C[i * ldc + j] = ddiconv(DIDEFAULTFOLD, CI + (i * N + j) * dinum(DIDEFAULTFOLD));
        }
      }
      break;
    default:
      if(beta == 1.0){
        for(j = 0; j < N; j++){
          for(i = 0; i < M; i++){
            didconv(DIDEFAULTFOLD, C[j * ldc + i], CI + (j * M + i) * dinum(DIDEFAULTFOLD));
          }
        }
      }else if(beta == 0.0){
        memset(CI, 0.0, M * N * disize(DIDEFAULTFOLD));
      }else{
        for(j = 0; j < N; j++){
          for(i = 0; i < M; i++){
            didconv(DIDEFAULTFOLD, C[j * ldc + i] * beta, CI + (j * M + i) * dinum(DIDEFAULTFOLD));
          }
        }
      }
      didgemm(DIDEFAULTFOLD, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, M);
      for(j = 0; j < N; j++){
        for(i = 0; i < M; i++){
          C[j * ldc + i] = ddiconv(DIDEFAULTFOLD, CI + (j * M + i) * dinum(DIDEFAULTFOLD));
        }
      }
      break;
  }
  free(CI);
}
