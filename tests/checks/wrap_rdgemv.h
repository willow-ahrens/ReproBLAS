#include <indexed.h>
#include <indexedBLAS.h>
#include <reproBLAS.h>
#include "../../config.h"

void wrap_rdgemv(int fold, const char Order, const char TransA, const int M, const int N, const double alpha, const double *A, const int lda, const double *X, const int incX, const double beta, double *Y, const int incY){
  double_indexed *YI;
  int NY;
  int i;
  if(fold == DIDEFAULTFOLD){
    rdgemv(Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
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
    YI = (double_indexed*)malloc(NY * incY * disize(fold));
    for(i = 0; i < NY; i++){
      didconv(fold, Y[i * incY] * beta, YI + i * incY * dinum(fold));
    }
    didgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, incY);
    for(i = 0; i < NY; i++){
      Y[i * incY] = ddiconv(fold, YI + i * incY * dinum(fold));
    }
    free(YI);
  }
}
