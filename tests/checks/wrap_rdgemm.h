#include <idxd.h>
#include <idxdBLAS.h>
#include <reproBLAS.h>

#include "../../config.h"

#include "wrap_daugsum.h"

void wrap_rdgemm(int fold, char Order, char TransA, char TransB, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc){
  double_indexed *CI;
  int i;
  int j;
  if(fold == DIDEFAULTFOLD){
    reproBLAS_rdgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }else{
    CI = (double_indexed*)malloc(M * N * idxd_disize(fold));
    if(beta == 0.0){
      memset(CI, 0, M * N * idxd_disize(fold));
    }else{
      for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
          switch(Order){
            case 'r':
            case 'R':
              idxd_didconv(fold, C[i * ldc + j] * beta, CI + (i * N + j) * idxd_dinum(fold));
              break;
            default:
              idxd_didconv(fold, C[i + j * ldc] * beta, CI + (i + j * M) * idxd_dinum(fold));
              break;
          }
        }
      }
    }
    idxdBLAS_didgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, ldc);
    for(i = 0; i < M; i++){
      for(j = 0; j < N; j++){
        switch(Order){
          case 'r':
          case 'R':
            C[i * ldc + j] = idxd_ddiconv(fold, CI + (i * N + j) * idxd_dinum(fold));
            break;
          default:
            C[i + j * ldc] = idxd_ddiconv(fold, CI + (i + j * M) * idxd_dinum(fold));
            break;
        }
      }
    }
    free(CI);
  }
}

void wrap_ref_rdgemm(int fold, char Order, char TransA, char TransB, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc){
  int i;
  int j;
  int k;
  double *opA = util_dmat_op(Order, TransA, M, K, A, lda);
  double *opB = util_dmat_op(Order, TransB, K, N, B, ldb);
  double_indexed *CI = idxd_dialloc(fold);
  if(alpha != 1.0){
    for(i = 0; i < M; i++){
      for(k = 0; k < K; k++){
        switch(Order){
          case 'r':
          case 'R':
            opA[i * K + k] = opA[i * K + k] * alpha;
            break;
          default:
            opA[k * M + i] = opA[k * M + i] * alpha;
            break;
        }
      }
    }
  }
  for(i = 0; i < M; i++){
    for(j = 0; j < N; j++){
      switch(Order){
        case 'r':
        case 'R':
          if(beta == 0.0){
            idxd_disetzero(fold, CI);
          }else{
            idxd_didconv(fold, C[i * ldc + j] * beta, CI);
          }
          if(alpha != 0.0){
            idxdBLAS_diddot(fold, K, opA + i * K, 1, opB + j, N, CI);
          }
          C[i * ldc + j] = idxd_ddiconv(fold, CI);
          break;
        default:
          if(beta == 0.0){
            idxd_disetzero(fold, CI);
          }else{
            idxd_didconv(fold, C[j * ldc + i] * beta, CI);
          }
          if(alpha != 0.0){
            idxdBLAS_diddot(fold, K, opA + i, M, opB + j * K, 1, CI);
          }
          C[j * ldc + i] = idxd_ddiconv(fold, CI);
          break;
      }
    }
  }
  free(CI);
  free(opA);
  free(opB);
}
