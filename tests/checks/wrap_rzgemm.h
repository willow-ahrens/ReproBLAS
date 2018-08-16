#include <binned.h>
#include <binnedBLAS.h>

#include "../common/test_util.h"

#include "../../config.h"

void wrap_rzgemm(int fold, char Order, char TransA, char TransB, int M, int N, int K, double complex *alpha, double complex *A, int lda, double complex *B, int ldb, double complex *beta, double complex *C, int ldc){
  if(fold == DIDEFAULTFOLD){
    reproBLAS_zgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }else{
    reproBLAS_rzgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  }
}

void wrap_ref_rzgemm(int fold, char Order, char TransA, char TransB, int M, int N, int K, double complex *alpha, double complex *A, int lda, double complex *B, int ldb, double complex *beta, double complex *C, int ldc){
  int i;
  int j;
  int k;
  double complex *opA = util_zmat_op(Order, TransA, M, K, A, lda);
  double complex *opB = util_zmat_op(Order, TransB, K, N, B, ldb);
  double complex betaC;
  double_complex_binned *CI = binned_zballoc(fold);
  if(*alpha != 1.0){
    for(i = 0; i < M; i++){
      for(k = 0; k < K; k++){
        switch(Order){
          case 'r':
          case 'R':
            opA[i * K + k] = zmul(opA[i * K + k], *alpha);
            break;
          default:
            opA[k * M + i] = zmul(opA[k * M + i], *alpha);
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
            binned_zbsetzero(fold, CI);
          }else if(*beta == 1.0){
            binned_zbzconv(fold, C + i * ldc + j, CI);
          }else{
            betaC = zmul(C[i * ldc + j], *beta);
            binned_zbzconv(fold, &betaC, CI);
          }
          if(*alpha != 0.0){
            binnedBLAS_zbzdotu(fold, K, opA + i * K, 1, opB + j, N, CI);
          }
          binned_zzbconv_sub(fold, CI, C + i * ldc + j);
          break;
        default:
          if(*beta == 0.0){
            binned_zbsetzero(fold, CI);
          }else if(*beta == 1.0){
            binned_zbzconv(fold, C + j * ldc + i, CI);
          }else{
            betaC = zmul(C[j * ldc + i], *beta);
            binned_zbzconv(fold, &betaC, CI);
          }
          if(*alpha != 0.0){
            binnedBLAS_zbzdotu(fold, K, opA + i, M, opB + j * K, 1, CI);
          }
          binned_zzbconv_sub(fold, CI, C + j * ldc + i);
          break;
      }
    }
  }
  free(CI);
  free(opA);
  free(opB);
}
