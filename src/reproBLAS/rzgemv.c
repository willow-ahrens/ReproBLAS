#include <stdlib.h>

#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

void rzgemv(const char Order,
            const char TransA, const int M, const int N,
            const void *alpha, const void *A, const int lda,
            const void *X, const int incX,
            const void *beta, void *Y, const int incY){
  double_complex_indexed *YI;
  double betaY[2];
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

  YI = (double_complex_indexed*)malloc(zisize(DIDEFAULTFOLD)*NY);
  for(i = 0; i < NY; i++){
    betaY[0] = ((double*)Y)[2 * i * incY] * ((double*)beta)[0] - ((double*)Y)[2 * i * incY + 1] * ((double*)beta)[1];
    betaY[1] = ((double*)Y)[2 * i * incY] * ((double*)beta)[1] + ((double*)Y)[2 * i * incY + 1] * ((double*)beta)[0];
    zizconv(DIDEFAULTFOLD, betaY, YI + i*zinum(DIDEFAULTFOLD));
  }
  zizgemv(DIDEFAULTFOLD, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
  for(i = 0; i < NY; i++){
    zziconv_sub(DIDEFAULTFOLD, YI + i*zinum(DIDEFAULTFOLD), ((double*)Y) + 2 * i * incY);
  }
  free(YI);
}
