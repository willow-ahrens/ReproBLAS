#include <stdlib.h>
#include <string.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

void reproBLAS_rdgemv(const int fold, const char Order,
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
      YI = (double_indexed*)malloc(M * idxd_disize(fold));
      if(beta == 0.0){
        memset(YI, 0, M * idxd_disize(fold));
      }else if(beta == 1.0){
        for(i = 0; i < M; i++){
          idxd_didconv(fold, Y[i * incY], YI + i * idxd_dinum(fold));
        }
      }else{
        for(i = 0; i < M; i++){
          idxd_didconv(fold, Y[i * incY] * beta, YI + i * idxd_dinum(fold));
        }
      }
      idxdBLAS_didgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < M; i++){
        Y[i * incY] = idxd_ddiconv(fold, YI + i * idxd_dinum(fold));
      }
      break;
    default:
      YI = (double_indexed*)malloc(N * idxd_disize(fold));
      if(beta == 0.0){
        memset(YI, 0, N * idxd_disize(fold));
      }else if(beta == 1.0){
        for(i = 0; i < N; i++){
          idxd_didconv(fold, Y[i * incY], YI + i * idxd_dinum(fold));
        }
      }else{
        for(i = 0; i < N; i++){
          idxd_didconv(fold, Y[i * incY] * beta, YI + i * idxd_dinum(fold));
        }
      }
      idxdBLAS_didgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < N; i++){
        Y[i * incY] = idxd_ddiconv(fold, YI + i * idxd_dinum(fold));
      }
      break;
  }
  free(YI);
}
