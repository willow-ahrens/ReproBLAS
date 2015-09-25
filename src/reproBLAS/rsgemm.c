#include <stdlib.h>
#include <string.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

void reproBLAS_rsgemm(const int fold, const char Order, const char TransA, const char TransB,
                      const int M, const int N, const int K,
                      const float alpha, const float *A, const int lda,
                      const float *B, const int ldb,
                      const float beta, float *C, const int ldc){
  float_indexed *CI;
  int i;
  int j;

  if(M == 0 || N == 0){
    return;
  }

  CI = (float_indexed*)malloc(M * N * idxd_sisize(fold));
  switch(Order){
    case 'r':
    case 'R':
      if(beta == 0.0){
        memset(CI, 0, M * N * idxd_sisize(fold));
      }else if(beta == 1.0){
        for(i = 0; i < M; i++){
          for(j = 0; j < N; j++){
            idxd_sisconv(fold, C[i * ldc + j], CI + (i * N + j) * idxd_sinum(fold));
          }
        }
      }else{
        for(i = 0; i < M; i++){
          for(j = 0; j < N; j++){
            idxd_sisconv(fold, C[i * ldc + j] * beta, CI + (i * N + j) * idxd_sinum(fold));
          }
        }
      }
      idxdBLAS_sisgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, N);
      for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
          C[i * ldc + j] = idxd_ssiconv(fold, CI + (i * N + j) * idxd_sinum(fold));
        }
      }
      break;
    default:
      if(beta == 0.0){
        memset(CI, 0, M * N * idxd_sisize(fold));
      }else if(beta == 1.0){
        for(j = 0; j < N; j++){
          for(i = 0; i < M; i++){
            idxd_sisconv(fold, C[j * ldc + i], CI + (j * M + i) * idxd_sinum(fold));
          }
        }
      }else{
        for(j = 0; j < N; j++){
          for(i = 0; i < M; i++){
            idxd_sisconv(fold, C[j * ldc + i] * beta, CI + (j * M + i) * idxd_sinum(fold));
          }
        }
      }
      idxdBLAS_sisgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, M);
      for(j = 0; j < N; j++){
        for(i = 0; i < M; i++){
          C[j * ldc + i] = idxd_ssiconv(fold, CI + (j * M + i) * idxd_sinum(fold));
        }
      }
      break;
  }
  free(CI);
}
