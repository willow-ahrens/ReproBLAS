#include <indexed.h>
#include <indexedBLAS.h>
#include <reproBLAS.h>
#include "../../config.h"

void wrap_rzgemv(int fold, char Order, char TransA, int M, int N, double complex *alpha, double complex *A, int lda, double complex *X, int incX, double complex *beta, double complex *Y, int incY){
  double_indexed *YI;
  double complex betaY;
  int NY;
  int i;
  if(fold == DIDEFAULTFOLD){
    rzgemv(Order, TransA, M, N, (void*)alpha, (void*)A, lda, (void*)X, incX, (void*)beta, (void*)Y, incY);
  }else{
    switch(TransA){
      case 'n':
      case 'N':
        NY = M;
      break;
      default:
        NY = N;
      break;
    }
    YI = (double_indexed*)malloc(NY * incY * zisize(fold));
    for(i = 0; i < NY; i++){
      betaY = Y[i * incY] * (*beta);
      zizconv(fold, &betaY, YI + i * incY * zinum(fold));
    }
    zizgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, incY);
    for(i = 0; i < NY; i++){
      zziconv_sub(fold, YI + i * incY * zinum(fold), Y + i * incY);
    }
    free(YI);
  }
}
