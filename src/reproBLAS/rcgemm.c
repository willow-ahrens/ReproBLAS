#include <stdlib.h>
#include <string.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

void reproBLAS_rcgemm(const int fold, const char Order, const char TransA, const char TransB,
                      const int M, const int N, const int K,
                      const void *alpha, const void *A, const int lda,
                      const void *B, const int ldb,
                      const void *beta, void *C, const int ldc){
  float_complex_indexed *CI;
  float betaC[2];
  int i;
  int j;

  if(M == 0 || N == 0){
    return;
  }

  CI = (float_complex_indexed*)malloc(M * N * idxd_cisize(fold));
  switch(Order){
    case 'r':
    case 'R':
      if(((float*)beta)[0] == 0.0 && ((float*)beta)[1] == 0.0){
        memset(CI, 0, M * N * idxd_cisize(fold));
      }else if(((float*)beta)[0] == 1.0 && ((float*)beta)[1] == 0.0){
        for(i = 0; i < M; i++){
          for(j = 0; j < N; j++){
            idxd_cicconv(fold, ((float*)C) + 2 * (i * ldc + j), CI + (i * N + j) * idxd_cinum(fold));
          }
        }
      }else{
        for(i = 0; i < M; i++){
          for(j = 0; j < N; j++){
            betaC[0] = ((float*)C)[2 * (i * ldc + j)] * ((float*)beta)[0] - ((float*)C)[2 * (i * ldc + j) + 1] * ((float*)beta)[1],
            betaC[1] = ((float*)C)[2 * (i * ldc + j)] * ((float*)beta)[1] + ((float*)C)[2 * (i * ldc + j) + 1] * ((float*)beta)[0],
            idxd_cicconv(fold, betaC, CI + (i * N + j) * idxd_cinum(fold));
          }
        }
      }
      idxdBLAS_cicgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, N);
      for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
          idxd_cciconv_sub(fold, CI + (i * N + j) * idxd_cinum(fold), ((float*)C) + 2 * (i * ldc + j));
        }
      }
      break;
    default:
      if(((float*)beta)[0] == 0.0 && ((float*)beta)[1] == 0.0){
        memset(CI, 0, M * N * idxd_cisize(fold));
      }else if(((float*)beta)[0] == 1.0 && ((float*)beta)[1] == 0.0){
        for(j = 0; j < N; j++){
          for(i = 0; i < M; i++){
            idxd_cicconv(fold, ((float*)C) + 2 * (j * ldc + i), CI + (j * M + i) * idxd_cinum(fold));
          }
        }
      }else{
        for(j = 0; j < N; j++){
          for(i = 0; i < M; i++){
            betaC[0] = ((float*)C)[2 * (j * ldc + i)] * ((float*)beta)[0] - ((float*)C)[2 * (j * ldc + i) + 1] * ((float*)beta)[1],
            betaC[1] = ((float*)C)[2 * (j * ldc + i)] * ((float*)beta)[1] + ((float*)C)[2 * (j * ldc + i) + 1] * ((float*)beta)[0],
            idxd_cicconv(fold, betaC, CI + (j * M + i) * idxd_cinum(fold));
          }
        }
      }
      idxdBLAS_cicgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, M);
      for(j = 0; j < N; j++){
        for(i = 0; i < M; i++){
          idxd_cciconv_sub(fold, CI + (j * M + i) * idxd_cinum(fold), ((float*)C) + 2 * (j * ldc + i));
        }
      }
      break;
  }
  free(CI);
}
