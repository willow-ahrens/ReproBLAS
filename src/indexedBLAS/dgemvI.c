#include "indexedBLAS.h"
#include "indexed.h"
#include "../Common/Common.h"

#define Y_BLOCK_SIZE 16
#define X_BLOCK_SIZE 1024

void dgemvI(rblas_order_t Order,
            rblas_transpose_t TransA, int M, int N,
            double *A, int lda,
            double *X, int incX,
            I_double *Y, int incY, int fold){
  switch(Order){
    case rblas_Row_Major:
      switch(TransA){
        case rblas_No_Trans:
          for(int i = 0; i < M; i += Y_BLOCK_SIZE){
            for(int j = 0; j < N; j += X_BLOCK_SIZE){
              for(int ii = i; ii < M && ii < i + Y_BLOCK_SIZE; ii++){
                ddotI1(MIN(X_BLOCK_SIZE, N - j), X + j, incX, A + ii * lda + j, 1, fold,&Y[ii*incY].m,&Y[ii*incY].c);
              }
            }
          }
          /*
          for(int i = 0; i < M; i++){
            ddotI1(N,X,incX, A + i * lda, 1, fold,&Y[i*incY].m,&Y[i*incY].c);
          }
          */
          break;
        default:
          for(int i = 0; i < N; i++){
            ddotI1(M,X,incX, A + i, lda, fold,&Y[i*incY].m,&Y[i*incY].c);
          }
          break;
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
