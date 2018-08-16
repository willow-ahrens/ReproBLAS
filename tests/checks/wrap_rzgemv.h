#include <binned.h>
#include <binnedBLAS.h>
#include <reproBLAS.h>

#include "../../config.h"

void wrap_rzgemv(int fold, char Order, char TransA, int M, int N, double complex *alpha, double complex *A, int lda, double complex *X, int incX, double complex *beta, double complex *Y, int incY){
  if(fold == DIDEFAULTFOLD){
    reproBLAS_zgemv(Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  }else{
    reproBLAS_rzgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  }
}

void wrap_ref_rzgemv(int fold, char Order, char TransA, int M, int N, double complex *alpha, double complex *A, int lda, double complex *X, int incX, double complex *beta, double complex *Y, int incY){
  int opM;
  int opN;
  double complex *opA;
  double complex *opX;
  double_complex_binned *YI;
  double complex betaY;
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
  opA = util_zmat_op(Order, TransA, opM, opN, A, lda);
  YI = binned_zballoc(fold);
  opX = (double complex*)malloc(opN * sizeof(double complex));
  for(i = 0; i < opM; i++){
    if(*beta == 0.0){
      binned_zbsetzero(fold, YI);
    }else if(*beta == 1.0){
      binned_zbzconv(fold, Y + i * incY, YI);
    }else{
      betaY = zmul(Y[i * incY], *beta);
      binned_zbzconv(fold, &betaY, YI);
    }
    if(*alpha != 0.0){
      for(j = 0; j < opN; j++){
        if(*alpha == 1.0){
          opX[j] = X[j * incX];
        }else{
          opX[j] = zmul(*alpha, X[j * incX]);
        }
      }
      switch(Order){
        case 'r':
        case 'R':
          binnedBLAS_zbzdotu(fold, opN, opA + i * opN, 1, opX, 1, YI);
          break;
        default:
          binnedBLAS_zbzdotu(fold, opN, opA + i, opM, opX, 1, YI);
          break;
      }
    }
    binned_zzbconv_sub(fold, YI, Y + i * incY);
  }
  free(YI);
  free(opA);
  free(opX);
}
