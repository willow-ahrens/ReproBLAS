#include <stdlib.h>
#include <string.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

void reproBLAS_rzgemm(const int fold, const char Order, const char TransA, const char TransB,
                      const int M, const int N, const int K,
                      const void *alpha, const void *A, const int lda,
                      const void *B, const int ldb,
                      const void *beta, void *C, const int ldc){
  double_complex_indexed *CI;
  double betaC[2];
  int i;
  int j;

  if(M == 0 || N == 0){
    return;
  }

  CI = (double_complex_indexed*)malloc(M * N * idxd_zisize(fold));
  switch(Order){
    case 'r':
    case 'R':
      if(((double*)beta)[0] == 0.0 && ((double*)beta)[1] == 0.0){
        memset(CI, 0, M * N * idxd_zisize(fold));
      }else if(((double*)beta)[0] == 1.0 && ((double*)beta)[1] == 0.0){
        for(i = 0; i < M; i++){
          for(j = 0; j < N; j++){
            idxd_zizconv(fold, ((double*)C) + 2 * (i * ldc + j), CI + (i * N + j) * idxd_zinum(fold));
          }
        }
      }else{
        for(i = 0; i < M; i++){
          for(j = 0; j < N; j++){
            betaC[0] = ((double*)C)[2 * (i * ldc + j)] * ((double*)beta)[0] - ((double*)C)[2 * (i * ldc + j) + 1] * ((double*)beta)[1],
            betaC[1] = ((double*)C)[2 * (i * ldc + j)] * ((double*)beta)[1] + ((double*)C)[2 * (i * ldc + j) + 1] * ((double*)beta)[0],
            idxd_zizconv(fold, betaC, CI + (i * N + j) * idxd_zinum(fold));
          }
        }
      }
      idxdBLAS_zizgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, N);
      for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
          idxd_zziconv_sub(fold, CI + (i * N + j) * idxd_zinum(fold), ((double*)C) + 2 * (i * ldc + j));
        }
      }
      break;
    default:
      if(((double*)beta)[0] == 0.0 && ((double*)beta)[1] == 0.0){
        memset(CI, 0, M * N * idxd_zisize(fold));
      }else if(((double*)beta)[0] == 1.0 && ((double*)beta)[1] == 0.0){
        for(j = 0; j < N; j++){
          for(i = 0; i < M; i++){
            idxd_zizconv(fold, ((double*)C) + 2 * (j * ldc + i), CI + (j * M + i) * idxd_zinum(fold));
          }
        }
      }else{
        for(j = 0; j < N; j++){
          for(i = 0; i < M; i++){
            betaC[0] = ((double*)C)[2 * (j * ldc + i)] * ((double*)beta)[0] - ((double*)C)[2 * (j * ldc + i) + 1] * ((double*)beta)[1],
            betaC[1] = ((double*)C)[2 * (j * ldc + i)] * ((double*)beta)[1] + ((double*)C)[2 * (j * ldc + i) + 1] * ((double*)beta)[0],
            idxd_zizconv(fold, betaC, CI + (j * M + i) * idxd_zinum(fold));
          }
        }
      }
      idxdBLAS_zizgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, M);
      for(j = 0; j < N; j++){
        for(i = 0; i < M; i++){
          idxd_zziconv_sub(fold, CI + (j * M + i) * idxd_zinum(fold), ((double*)C) + 2 * (j * ldc + i));
        }
      }
      break;
  }
  free(CI);
}
