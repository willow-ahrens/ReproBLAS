#include "indexedBLAS.h"
#include "indexed.h"

void dgemvI(rblas_order_t Order,
            rblas_transpose_t TransA, int M, int N,
            double *A, int lda,
            double *X, int incX,
            I_double *Y, int incY, int fold){
  switch(Order){
    case rblas_Row_Major:
      switch(TransA){
        case rblas_No_Trans:
          for(int i = 0; i < M; i++){
            ddotI1(N,X,incX, A + i * lda, 1, fold,40,&Y[i*incY].m,&Y[i*incY].c);
          }
          break;
        default:
          for(int i = 0; i < N; i++){
            ddotI1(M,X,incX, A + i, lda, fold,40,&Y[i*incY].m,&Y[i*incY].c);
          }
      }
      break;
    case rblas_Col_Major:
      switch(TransA){
        case rblas_No_Trans:
          dgemvI(rblas_Row_Major, rblas_Trans, N, M, A, lda, X, incX, Y, incY, fold);
          break;
        default:
          dgemvI(rblas_Row_Major, rblas_No_Trans, N, M, A, lda, X, incX, Y, incY, fold);
          break;
      }
      break;
  }
}
