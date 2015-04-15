#include <indexedBLAS.h>
#include <reproBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_util.h"

int main(int argc, char**argv){
  util_random_seed();
  int N = 4;
  int M = 5;
  int lda = 6;
  double *X = util_dvec_alloc(N, 1);
  double *Y = util_dvec_alloc(M, 1);
  double *A = util_dmat_alloc(rblas_Row_Major, M, N,lda);
  util_dvec_fill(N, X, 1, util_Vec_Rand, 1, 1);
  util_dvec_fill(M, Y, 1, util_Vec_Constant, 0, 1);
  util_dmat_fill(rblas_Row_Major, rblas_No_Trans, M, N, A, lda, util_Mat_Row_Rand, 1, 1);

  rdgemv(rblas_Row_Major,
         rblas_No_Trans, M, N,
         A, lda,
         X, 1,
         Y, 1);

  printf("Y X:\t");
  for(int i = 0; i < N; i++){
    printf("%.5f\t", X[i]);
  }
  printf("\n");
  for(int i = 0; i < M; i++){
    printf("%.5f\t", Y[i]);
    for(int j = 0; j < N; j++){
      printf("%.5f\t", A[i * lda + j]);
    }
    printf("\n");
  }
  free(X);
  free(Y);
  free(A);
}
