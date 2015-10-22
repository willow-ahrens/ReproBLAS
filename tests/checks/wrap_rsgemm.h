#include <idxd.h>
#include <idxdBLAS.h>

#include "../common/test_util.h"

#include "../../config.h"

void wrap_rsgemm(int fold, char Order, char TransA, char TransB, int M, int N, int K, float alpha, float *A, int lda, float *B, int ldb, float beta, float *C, int ldc){
  if(fold == SIDEFAULTFOLD){
    reproBLAS_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }else{
    reproBLAS_rsgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }
}

void wrap_ref_rsgemm(int fold, char Order, char TransA, char TransB, int M, int N, int K, float alpha, float *A, int lda, float *B, int ldb, float beta, float *C, int ldc){
  int i;
  int j;
  int k;
  float *opA = util_smat_op(Order, TransA, M, K, A, lda);
  float *opB = util_smat_op(Order, TransB, K, N, B, ldb);
  float_indexed *CI = idxd_sialloc(fold);
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
            idxd_sisetzero(fold, CI);
          }else{
            idxd_sisconv(fold, C[i * ldc + j] * beta, CI);
          }
          if(alpha != 0.0){
            idxdBLAS_sisdot(fold, K, opA + i * K, 1, opB + j, N, CI);
          }
          C[i * ldc + j] = idxd_ssiconv(fold, CI);
          break;
        default:
          if(beta == 0.0){
            idxd_sisetzero(fold, CI);
          }else{
            idxd_sisconv(fold, C[j * ldc + i] * beta, CI);
          }
          if(alpha != 0.0){
            idxdBLAS_sisdot(fold, K, opA + i, M, opB + j * K, 1, CI);
          }
          C[j * ldc + i] = idxd_ssiconv(fold, CI);
          break;
      }
    }
  }
  free(CI);
  free(opA);
  free(opB);
}
