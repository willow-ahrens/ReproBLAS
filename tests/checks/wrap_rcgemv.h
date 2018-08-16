#include <binned.h>
#include <binnedBLAS.h>
#include <reproBLAS.h>

#include "../../config.h"

void wrap_rcgemv(int fold, char Order, char TransA, int M, int N, float complex *alpha, float complex *A, int lda, float complex *X, int incX, float complex *beta, float complex *Y, int incY){
  if(fold == SIDEFAULTFOLD){
    reproBLAS_cgemv(Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  }else{
    reproBLAS_rcgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  }
}

void wrap_ref_rcgemv(int fold, char Order, char TransA, int M, int N, float complex *alpha, float complex *A, int lda, float complex *X, int incX, float complex *beta, float complex *Y, int incY){
  int opM;
  int opN;
  float complex *opA;
  float complex *opX;
  float_complex_binned *YI;
  float complex betaY;
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
  opA = util_cmat_op(Order, TransA, opM, opN, A, lda);
  YI = binned_cballoc(fold);
  opX = (float complex*)malloc(opN * sizeof(float complex));
  for(i = 0; i < opM; i++){
    if(*beta == 0.0){
      binned_cbsetzero(fold, YI);
    }else if(*beta == 1.0){
      binned_cbcconv(fold, Y + i * incY, YI);
    }else{
      betaY = cmul(Y[i * incY], *beta);
      binned_cbcconv(fold, &betaY, YI);
    }
    if(*alpha != 0.0){
      for(j = 0; j < opN; j++){
        if(*alpha == 1.0){
          opX[j] = X[j * incX];
        }else{
          opX[j] = cmul(*alpha, X[j * incX]);
        }
      }
      switch(Order){
        case 'r':
        case 'R':
          binnedBLAS_cbcdotu(fold, opN, opA + i * opN, 1, opX, 1, YI);
          break;
        default:
          binnedBLAS_cbcdotu(fold, opN, opA + i, opM, opX, 1, YI);
          break;
      }
    }
    binned_ccbconv_sub(fold, YI, Y + i * incY);
  }
  free(YI);
  free(opA);
  free(opX);
}
