#include <idxd.h>
#include <idxdBLAS.h>

#include "../common/test_util.h"

#include "../../config.h"

void wrap_rdgemm(int fold, char Order, char TransA, char TransB, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc){
  if(fold == DIDEFAULTFOLD){
    reproBLAS_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }else{
    reproBLAS_rdgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
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
