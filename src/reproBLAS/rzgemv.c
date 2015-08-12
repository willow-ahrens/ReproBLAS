#include <stdlib.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

#include "../../config.h"

void rzgemv(const char Order,
            const char TransA, const int M, const int N,
            const void *alpha, const void *A, const int lda,
            const void *X, const int incX,
            const void *beta, void *Y, const int incY){
  double_complex_indexed *YI;
  double betaY[2];
  int i;

  switch(TransA){
    case 'n':
    case 'N':
      YI = (double_complex_indexed*)malloc(idxd_zisize(DIDEFAULTFOLD)*M);
      for(i = 0; i < M; i++){
        betaY[0] = ((double*)Y)[2 * i * incY] * ((double*)beta)[0] - ((double*)Y)[2 * i * incY + 1] * ((double*)beta)[1];
        betaY[1] = ((double*)Y)[2 * i * incY] * ((double*)beta)[1] + ((double*)Y)[2 * i * incY + 1] * ((double*)beta)[0];
        idxd_zizconv(DIDEFAULTFOLD, betaY, YI + i * idxd_zinum(DIDEFAULTFOLD));
      }
      idxdBLAS_zizgemv(DIDEFAULTFOLD, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < M; i++){
        idxd_zziconv_sub(DIDEFAULTFOLD, YI + i * idxd_zinum(DIDEFAULTFOLD), ((double*)Y) + 2 * i * incY);
      }
      break;
    default:
      YI = (double_complex_indexed*)malloc(idxd_zisize(DIDEFAULTFOLD)*N);
      for(i = 0; i < N; i++){
        betaY[0] = ((double*)Y)[2 * i * incY] * ((double*)beta)[0] - ((double*)Y)[2 * i * incY + 1] * ((double*)beta)[1];
        betaY[1] = ((double*)Y)[2 * i * incY] * ((double*)beta)[1] + ((double*)Y)[2 * i * incY + 1] * ((double*)beta)[0];
        idxd_zizconv(DIDEFAULTFOLD, betaY, YI + i * idxd_zinum(DIDEFAULTFOLD));
      }
      idxdBLAS_zizgemv(DIDEFAULTFOLD, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < N; i++){
        idxd_zziconv_sub(DIDEFAULTFOLD, YI + i * idxd_zinum(DIDEFAULTFOLD), ((double*)Y) + 2 * i * incY);
      }
      break;
  }

  free(YI);
}
