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
          YI = (Idouble*)malloc(dISize(DEFAULT_FOLD)*M);
          for(int i = 0; i < M; i++){
            YI[i] = dconv2I(Y[i * incY]);
          }
          dgemvI(order, TransA, M, N, A, lda, X, incX, YI, 1, DEFAULT_FOLD);
          for(int i = 0; i < M; i++){
            Y[i * incY] = Iconv2d1(DEFAULT_FOLD, &YI[i].m, &YI[i].c, 1);
          }
          free(YI);
          break;
        default:
          YI = (Idouble*)malloc(dISize(DEFAULT_FOLD)*N);
          for(int i = 0; i < N; i++){
            YI[i] = dconv2I(Y[i * incY]);
          }
          dgemvI(order, TransA, M, N, A, lda, X, incX, YI, 1, DEFAULT_FOLD);
          for(int i = 0; i < N; i++){
            Y[i * incY] = Iconv2d1(DEFAULT_FOLD, &YI[i].m, &YI[i].c, 1);
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
