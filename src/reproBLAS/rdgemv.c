#include <stdlib.h>

#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

void rdgemv(const char Order,
            const char TransA, const int M, const int N,
            const double alpha, const double *A, const int lda,
            const double *X, const int incX,
            const double beta, double *Y, const int incY){
  double_indexed *YI;
  int NY;
  int i;

  switch(TransA){
    case 'n':
    case 'N':
      NY = M;
      break;
    default:
      NY = N;
      break;
  }

  YI = (double_indexed*)malloc(disize(DIDEFAULTFOLD)*NY);
  for(i = 0; i < NY; i++){
    didconv(DIDEFAULTFOLD, Y[i * incY] * beta, YI + i*dinum(DIDEFAULTFOLD));
  }
  didgemv(DIDEFAULTFOLD, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
  for(i = 0; i < NY; i++){
    Y[i * incY] = ddiconv(DIDEFAULTFOLD, YI + i*dinum(DIDEFAULTFOLD));
  }
  free(YI);
}
