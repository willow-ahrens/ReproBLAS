#include <stdlib.h>
#include <string.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

void reproBLAS_rcgemv(const int fold, const char Order,
                      const char TransA, const int M, const int N,
                      const void *alpha, const void *A, const int lda,
                      const void *X, const int incX,
                      const void *beta, void *Y, const int incY){
  float_complex_indexed *YI;
  float betaY[2];
  int i;

  if(N == 0 || M == 0){
    return;
  }

  switch(TransA){
    case 'n':
    case 'N':
      YI = (float_complex_indexed*)malloc(M * idxd_cisize(fold));
      if(((float*)beta)[0] == 0.0 && ((float*)beta)[1] == 0.0){
        memset(YI, 0, M * idxd_cisize(fold));
      }else if(((float*)beta)[0] == 1.0 && ((float*)beta)[1] == 0.0){
        for(i = 0; i < M; i++){
          idxd_cicconv(fold, ((float*)Y) + 2 * i * incY, YI + i * idxd_cinum(fold));
        }
      }else{
        for(i = 0; i < M; i++){
          betaY[0] = ((float*)Y)[2 * i * incY] * ((float*)beta)[0] - ((float*)Y)[2 * i * incY + 1] * ((float*)beta)[1];
          betaY[1] = ((float*)Y)[2 * i * incY] * ((float*)beta)[1] + ((float*)Y)[2 * i * incY + 1] * ((float*)beta)[0];
          idxd_cicconv(fold, betaY, YI + i * idxd_cinum(fold));
        }
      }
      idxdBLAS_cicgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < M; i++){
        idxd_cciconv_sub(fold, YI + i * idxd_cinum(fold), ((float*)Y) + 2 * i * incY);
      }
      break;
    default:
      YI = (float_complex_indexed*)malloc(N * idxd_cisize(fold));
      if(((float*)beta)[0] == 0.0 && ((float*)beta)[1] == 0.0){
        memset(YI, 0, N * idxd_cisize(fold));
      }else if(((float*)beta)[0] == 1.0 && ((float*)beta)[1] == 0.0){
        for(i = 0; i < N; i++){
          idxd_cicconv(fold, ((float*)Y) + 2 * i * incY, YI + i * idxd_cinum(fold));
        }
      }else{
        for(i = 0; i < N; i++){
          betaY[0] = ((float*)Y)[2 * i * incY] * ((float*)beta)[0] - ((float*)Y)[2 * i * incY + 1] * ((float*)beta)[1];
          betaY[1] = ((float*)Y)[2 * i * incY] * ((float*)beta)[1] + ((float*)Y)[2 * i * incY + 1] * ((float*)beta)[0];
          idxd_cicconv(fold, betaY, YI + i * idxd_cinum(fold));
        }
      }
      idxdBLAS_cicgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < N; i++){
        idxd_cciconv_sub(fold, YI + i * idxd_cinum(fold), ((float*)Y) + 2 * i * incY);
      }
      break;
  }

  free(YI);
}
