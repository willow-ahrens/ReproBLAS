#include <stdlib.h>
#include <string.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

void reproBLAS_rsgemv(const int fold, const char Order,
                      const char TransA, const int M, const int N,
                      const float alpha, const float *A, const int lda,
                      const float *X, const int incX,
                      const float beta, float *Y, const int incY){
  float_indexed *YI;
  int i;

  if(N == 0 || M == 0){
    return;
  }

  switch(TransA){
    case 'n':
    case 'N':
      YI = (float_indexed*)malloc(M * idxd_sisize(fold));
      if(beta == 0.0){
        memset(YI, 0, M * idxd_sisize(fold));
      }else if(beta == 1.0){
        for(i = 0; i < M; i++){
          idxd_sisconv(fold, Y[i * incY], YI + i * idxd_sinum(fold));
        }
      }else{
        for(i = 0; i < M; i++){
          idxd_sisconv(fold, Y[i * incY] * beta, YI + i * idxd_sinum(fold));
        }
      }
      idxdBLAS_sisgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < M; i++){
        Y[i * incY] = idxd_ssiconv(fold, YI + i * idxd_sinum(fold));
      }
      break;
    default:
      YI = (float_indexed*)malloc(N * idxd_sisize(fold));
      if(beta == 0.0){
        memset(YI, 0, N * idxd_sisize(fold));
      }else if(beta == 1.0){
        for(i = 0; i < N; i++){
          idxd_sisconv(fold, Y[i * incY], YI + i * idxd_sinum(fold));
        }
      }else{
        for(i = 0; i < N; i++){
          idxd_sisconv(fold, Y[i * incY] * beta, YI + i * idxd_sinum(fold));
        }
      }
      idxdBLAS_sisgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < N; i++){
        Y[i * incY] = idxd_ssiconv(fold, YI + i * idxd_sinum(fold));
      }
      break;
  }
  free(YI);
}
