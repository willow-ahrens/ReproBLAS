#include <stdlib.h>
#include <string.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

void reproBLAS_rzgemv(const int fold, const char Order,
            const char TransA, const int M, const int N,
            const void *alpha, const void *A, const int lda,
            const void *X, const int incX,
            const void *beta, void *Y, const int incY){
  double_complex_indexed *YI;
  double betaY[2];
  int i;

  if(N == 0 || M == 0){
    return;
  }

  switch(TransA){
    case 'n':
    case 'N':
      YI = (double_complex_indexed*)malloc(M * idxd_zisize(fold));
      if(((double*)beta)[0] == 0.0 && ((double*)beta)[1] == 0.0){
        memset(YI, 0, M * idxd_zisize(fold));
      }else if(((double*)beta)[0] == 1.0 && ((double*)beta)[1] == 0.0){
        for(i = 0; i < M; i++){
          idxd_zizconv(fold, ((double*)Y) + 2 * i * incY, YI + i * idxd_zinum(fold));
        }
      }else{
        for(i = 0; i < M; i++){
          betaY[0] = ((double*)Y)[2 * i * incY] * ((double*)beta)[0] - ((double*)Y)[2 * i * incY + 1] * ((double*)beta)[1];
          betaY[1] = ((double*)Y)[2 * i * incY] * ((double*)beta)[1] + ((double*)Y)[2 * i * incY + 1] * ((double*)beta)[0];
          idxd_zizconv(fold, betaY, YI + i * idxd_zinum(fold));
        }
      }
      idxdBLAS_zizgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < M; i++){
        idxd_zziconv_sub(fold, YI + i * idxd_zinum(fold), ((double*)Y) + 2 * i * incY);
      }
      break;
    default:
      YI = (double_complex_indexed*)malloc(N * idxd_zisize(fold));
      if(((double*)beta)[0] == 0.0 && ((double*)beta)[1] == 0.0){
        memset(YI, 0, N * idxd_zisize(fold));
      }else if(((double*)beta)[0] == 1.0 && ((double*)beta)[1] == 0.0){
        for(i = 0; i < N; i++){
          idxd_zizconv(fold, ((double*)Y) + 2 * i * incY, YI + i * idxd_zinum(fold));
        }
      }else{
        for(i = 0; i < N; i++){
          betaY[0] = ((double*)Y)[2 * i * incY] * ((double*)beta)[0] - ((double*)Y)[2 * i * incY + 1] * ((double*)beta)[1];
          betaY[1] = ((double*)Y)[2 * i * incY] * ((double*)beta)[1] + ((double*)Y)[2 * i * incY + 1] * ((double*)beta)[0];
          idxd_zizconv(fold, betaY, YI + i * idxd_zinum(fold));
        }
      }
      idxdBLAS_zizgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < N; i++){
        idxd_zziconv_sub(fold, YI + i * idxd_zinum(fold), ((double*)Y) + 2 * i * incY);
      }
      break;
  }

  free(YI);
}
