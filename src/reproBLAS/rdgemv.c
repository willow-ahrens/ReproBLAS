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
  int i;

  if(N == 0 || M == 0){
    return;
  }

  switch(TransA){
    case 'n':
    case 'N':
      YI = (double_indexed*)malloc(M * idxd_disize(DIDEFAULTFOLD));
      if(beta == 1.0){
        for(i = 0; i < M; i++){
          idxd_didconv(DIDEFAULTFOLD, Y[i * incY], YI + i * idxd_dinum(DIDEFAULTFOLD));
        }
      }else if(beta == 0.0){
        memset(YI, 0, M * idxd_disize(DIDEFAULTFOLD));
      }else{
        for(i = 0; i < M; i++){
          idxd_didconv(DIDEFAULTFOLD, Y[i * incY] * beta, YI + i * idxd_dinum(DIDEFAULTFOLD));
        }
      }
      idxdBLAS_didgemv(DIDEFAULTFOLD, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < M; i++){
        Y[i * incY] = idxd_ddiconv(DIDEFAULTFOLD, YI + i * idxd_dinum(DIDEFAULTFOLD));
      }
      break;
    default:
      YI = (double_indexed*)malloc(N * idxd_disize(DIDEFAULTFOLD));
      if(beta == 1.0){
        for(i = 0; i < N; i++){
          idxd_didconv(DIDEFAULTFOLD, Y[i * incY], YI + i * idxd_dinum(DIDEFAULTFOLD));
        }
      }else if(beta == 0.0){
        memset(YI, 0, N * idxd_disize(DIDEFAULTFOLD));
      }else{
        for(i = 0; i < N; i++){
          idxd_didconv(DIDEFAULTFOLD, Y[i * incY] * beta, YI + i * idxd_dinum(DIDEFAULTFOLD));
        }
      }
      idxdBLAS_didgemv(DIDEFAULTFOLD, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < N; i++){
        Y[i * incY] = idxd_ddiconv(DIDEFAULTFOLD, YI + i * idxd_dinum(DIDEFAULTFOLD));
      }
      break;
  }
  free(YI);
}
