#include "indexedBLAS.h"
#include "indexed.h"
#include "../common/common.h"

#define Y_BLOCK_SIZE 16
#define X_BLOCK_SIZE 1024

void dgemvI(const int fold, const rblas_order_t Order,
            const rblas_transpose_t TransA, const int M, const int N,
            const double *A, const int lda,
            const double *X, const int incX,
            double_indexed *Y, const int incY){
  int i;
  int ii;
  int j;
  switch(Order){
    case rblas_Row_Major:
      switch(TransA){
        case rblas_No_Trans:
          for(i = 0; i < M; i += Y_BLOCK_SIZE){
            for(j = 0; j < N; j += X_BLOCK_SIZE){
              for(ii = i; ii < M && ii < i + Y_BLOCK_SIZE; ii++){
                diddot(fold, MIN(X_BLOCK_SIZE, N - j), X + j, incX, A + ii * lda + j, 1, Y + ii*incY*dinum(fold));
              }
            }
          }
          break;
        default:
          for(i = 0; i < N; i++){
            diddot(fold, M,X,incX, A + i, lda, Y + i*incY*dinum(fold));
          }
          break;
      }
      break;
    case rblas_Col_Major:
      switch(TransA){
        case rblas_No_Trans:
          dgemvI(fold, rblas_Row_Major, rblas_Trans, N, M, A, lda, X, incX, Y, incY);
          break;
        default:
          dgemvI(fold, rblas_Row_Major, rblas_No_Trans, N, M, A, lda, X, incX, Y, incY);
          break;
      }
      break;
  }
}
