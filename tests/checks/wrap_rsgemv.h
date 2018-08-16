#include <binned.h>
#include <binnedBLAS.h>

#include "../common/test_util.h"

#include "../../config.h"

void wrap_rsgemv(int fold, char Order, char TransA, int M, int N, float alpha, float *A, int lda, float *X, int incX, float beta, float *Y, int incY){
  if(fold == SIDEFAULTFOLD){
    reproBLAS_sgemv(Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  }else{
    reproBLAS_rsgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  }
}

void wrap_ref_rsgemv(int fold, char Order, char TransA, int M, int N, float alpha, float *A, int lda, float *X, int incX, float beta, float *Y, int incY){
  int opM;
  int opN;
  float *opA;
  float *opX;
  float_binned *YI;
  int i;
  int j;
  switch(TransA){
    case 'n':
    case 'N':
      opM = M;
      opN = N;
      break;
    default:
      opM = N;
      opN = M;
      break;
  }
  opA = util_smat_op(Order, TransA, opM, opN, A, lda);
  YI = binned_sballoc(fold);
  opX = (float*)malloc(opN * sizeof(float));
  for(i = 0; i < opM; i++){
    if(beta == 0.0){
      binned_sbsetzero(fold, YI);
    }else{
      binned_sbsconv(fold, Y[i * incY] * beta, YI);
    }
    if(alpha != 0.0){
      for(j = 0; j < opN; j++){
        opX[j] = alpha * X[j * incX];
      }
      switch(Order){
        case 'r':
        case 'R':
          binnedBLAS_sbsdot(fold, opN, opA + i * opN, 1, opX, 1, YI);
          break;
        default:
          binnedBLAS_sbsdot(fold, opN, opA + i, opM, opX, 1, YI);
          break;
      }
    }
    Y[i * incY] = binned_ssbconv(fold, YI);
  }
  free(YI);
  free(opA);
  free(opX);
}
