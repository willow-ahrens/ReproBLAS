#include <stdlib.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

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

  CI = (double_indexed*)malloc(M * N * idxd_disize(DIDEFAULTFOLD));
  switch(Order){
    case 'r':
    case 'R':
      if(beta == 1.0){
        for(i = 0; i < M; i++){
          for(j = 0; j < N; j++){
            idxd_didconv(DIDEFAULTFOLD, C[i * ldc + j], CI + (i * N + j) * idxd_dinum(DIDEFAULTFOLD));
          }
        }
      }else if(beta == 0.0){
        memset(CI, 0.0, M * N * idxd_disize(DIDEFAULTFOLD));
      }else{
        for(i = 0; i < M; i++){
          for(j = 0; j < N; j++){
            idxd_didconv(DIDEFAULTFOLD, C[i * ldc + j] * beta, CI + (i * N + j) * idxd_dinum(DIDEFAULTFOLD));
          }
        }
      }
      idxdBLAS_didgemm(DIDEFAULTFOLD, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, N);
      for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
          C[i * ldc + j] = idxd_ddiconv(DIDEFAULTFOLD, CI + (i * N + j) * idxd_dinum(DIDEFAULTFOLD));
        }
      }
      break;
    default:
      if(beta == 1.0){
        for(j = 0; j < N; j++){
          for(i = 0; i < M; i++){
            idxd_didconv(DIDEFAULTFOLD, C[j * ldc + i], CI + (j * M + i) * idxd_dinum(DIDEFAULTFOLD));
          }
        }
      }else if(beta == 0.0){
        memset(CI, 0.0, M * N * idxd_disize(DIDEFAULTFOLD));
      }else{
        for(j = 0; j < N; j++){
          for(i = 0; i < M; i++){
            idxd_didconv(DIDEFAULTFOLD, C[j * ldc + i] * beta, CI + (j * M + i) * idxd_dinum(DIDEFAULTFOLD));
          }
        }
      }
      idxdBLAS_didgemm(DIDEFAULTFOLD, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, M);
      for(j = 0; j < N; j++){
        for(i = 0; i < M; i++){
          C[j * ldc + i] = idxd_ddiconv(DIDEFAULTFOLD, CI + (j * M + i) * idxd_dinum(DIDEFAULTFOLD));
        }
      }
      break;
  }
  free(CI);
}
