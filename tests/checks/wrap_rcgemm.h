#include <idxd.h>
#include <idxdBLAS.h>

#include "../common/test_util.h"

#include "../../config.h"

void wrap_rcgemm(int fold, char Order, char TransA, char TransB, int M, int N, int K, float complex *alpha, float complex *A, int lda, float complex *B, int ldb, float complex *beta, float complex *C, int ldc){
  if(fold == SIDEFAULTFOLD){
    reproBLAS_cgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }else{
    reproBLAS_rcgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }
}

void wrap_ref_rcgemm(int fold, char Order, char TransA, char TransB, int M, int N, int K, float complex *alpha, float complex *A, int lda, float complex *B, int ldb, float complex *beta, float complex *C, int ldc){
  int i;
  int j;
  int k;
  float complex *opA = util_cmat_op(Order, TransA, M, K, A, lda);
  float complex *opB = util_cmat_op(Order, TransB, K, N, B, ldb);
  float complex betaC;
  float_complex_indexed *CI = idxd_cialloc(fold);
  if(*alpha != 1.0){
    for(i = 0; i < M; i++){
      for(k = 0; k < K; k++){
        switch(Order){
          case 'r':
          case 'R':
            opA[i * K + k] = cmul(opA[i * K + k], *alpha);
            break;
          default:
            opA[k * M + i] = cmul(opA[k * M + i], *alpha);
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
          if(*beta == 0.0){
            idxd_cisetzero(fold, CI);
          }else if(*beta == 1.0){
            idxd_cicconv(fold, C + i * ldc + j, CI);
          }else{
            betaC = cmul(C[i * ldc + j], *beta);
            idxd_cicconv(fold, &betaC, CI);
          }
          if(*alpha != 0.0){
            idxdBLAS_cicdotu(fold, K, opA + i * K, 1, opB + j, N, CI);
          }
          idxd_cciconv_sub(fold, CI, C + i * ldc + j);
          break;
        default:
          if(*beta == 0.0){
            idxd_cisetzero(fold, CI);
          }else if(*beta == 1.0){
            idxd_cicconv(fold, C + j * ldc + i, CI);
          }else{
            betaC = cmul(C[j * ldc + i], *beta);
            idxd_cicconv(fold, &betaC, CI);
          }
          if(*alpha != 0.0){
            idxdBLAS_cicdotu(fold, K, opA + i, M, opB + j * K, 1, CI);
          }
          idxd_cciconv_sub(fold, CI, C + j * ldc + i);
          break;
      }
    }
  }
  free(CI);
  free(opA);
  free(opB);
}
