/*
 *  Created   February 2015 Peter Ahrens
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"
#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

void rdgemv(const rblas_order_t order,
            const rblas_transpose_t TransA, const int M, const int N,
            const double *A, const int lda,
            const double *X, const int incX,
            double *Y, const int incY){
  double_indexed *YI;
  int i;

  switch(order){
    case rblas_Row_Major:
      switch(TransA){
        case rblas_No_Trans:
          YI = (double_indexed*)malloc(disize(DIDEFAULTFOLD)*M);
          for(i = 0; i < M; i++){
            didconv(DIDEFAULTFOLD, Y[i * incY], YI + i*dinum(DIDEFAULTFOLD));
          }
          dgemvI(DIDEFAULTFOLD, order, TransA, M, N, A, lda, X, incX, YI, 1);
          for(i = 0; i < M; i++){
            Y[i * incY] = ddiconv(DIDEFAULTFOLD, YI + i*dinum(DIDEFAULTFOLD));
          }
          free(YI);
          break;
        default:
          YI = (double_indexed*)malloc(disize(DIDEFAULTFOLD)*N);
          for(i = 0; i < N; i++){
            didconv(DIDEFAULTFOLD, Y[i * incY], YI + i*dinum(DIDEFAULTFOLD));
          }
          dgemvI(DIDEFAULTFOLD, order, TransA, M, N, A, lda, X, incX, YI, 1);
          for(i = 0; i < N; i++){
            Y[i * incY] = ddiconv(DIDEFAULTFOLD, YI + i*dinum(DIDEFAULTFOLD));
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
