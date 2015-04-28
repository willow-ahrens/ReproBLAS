/*
 *  Created   February 2015 Peter Ahrens
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

void rdgemv(const rblas_order_t order,
            const rblas_transpose_t TransA, const int M, const int N,
            const double *A, const int lda,
            const double *X, const int incX,
            double *Y, const int incY){
  Idouble *YI;
  switch(order){
    case rblas_Row_Major:
      switch(TransA){
        case rblas_No_Trans:
          YI = (Idouble*)malloc(disize(DEFAULT_FOLD)*M);
          for(int i = 0; i < M; i++){
            didconv(Y[i * incY], YI + i, DEFAULT_FOLD);
          }
          dgemvI(order, TransA, M, N, A, lda, X, incX, YI, 1, DEFAULT_FOLD);
          for(int i = 0; i < M; i++){
            Y[i * incY] = ddiconv(YI + i, DEFAULT_FOLD);
          }
          free(YI);
          break;
        default:
          YI = (Idouble*)malloc(disize(DEFAULT_FOLD)*N);
          for(int i = 0; i < N; i++){
            didconv(Y[i * incY], YI + i, DEFAULT_FOLD);
          }
          dgemvI(order, TransA, M, N, A, lda, X, incX, YI, 1, DEFAULT_FOLD);
          for(int i = 0; i < N; i++){
            Y[i * incY] = ddiconv(YI + i, DEFAULT_FOLD);
          }
          free(YI);
          break;
      }
      break;
    case rblas_Col_Major:
      switch(TransA){
        case rblas_No_Trans:
          rdgemv(rblas_Row_Major, rblas_Trans, N, M, A, lda, X, incX, Y, incY);
          break;
        default:
          rdgemv(rblas_Row_Major, rblas_No_Trans, N, M, A, lda, X, incX, Y, incY);
          break;
      }
      break;
  }
}
