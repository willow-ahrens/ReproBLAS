#include <stdlib.h>
#include <string.h>

#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Add to single precision vector Y the reproducible matrix-vector product of single precision matrix A and single precision vector X
 *
 * Performs one of the matrix-vector operations
 *
 *   y := alpha*A*x + beta*y   or   y := alpha*A**T*x + beta*y,
 *
 * where alpha and beta are scalars, x and y are vectors, and A is an M by N matrix.
 *
 * The matrix-vector product is computed using binned types with #binnedBLAS_sbsgemv()
 *
 * @param fold the fold of the binned types
 * @param Order a character specifying the matrix ordering ('r' or 'R' for row-major, 'c' or 'C' for column major)
 * @param TransA a character specifying whether or not to transpose A before taking the matrix-vector product ('n' or 'N' not to transpose, 't' or 'T' or 'c' or 'C' to transpose)
 * @param M number of rows of matrix A
 * @param N number of columns of matrix A
 * @param alpha scalar alpha
 * @param A single precision matrix of dimension (M, lda) in row-major or (lda, N) in column-major
 * @param lda the first dimension of A as declared in the calling program
 * @param X single precision vector of at least size N if not transposed or size M otherwise
 * @param incX X vector stride (use every incX'th element)
 * @param beta scalar beta
 * @param Y single precision vector Y of at least size M if not transposed or size N otherwise
 * @param incY Y vector stride (use every incY'th element)
 *
 * @author Willow Ahrens
 * @date   18 Jan 2016
 */
void reproBLAS_rsgemv(const int fold, const char Order,
                      const char TransA, const int M, const int N,
                      const float alpha, const float *A, const int lda,
                      const float *X, const int incX,
                      const float beta, float *Y, const int incY){
  float_binned *YI;
  int i;

  if(N == 0 || M == 0){
    return;
  }

  switch(TransA){
    case 'n':
    case 'N':
      YI = (float_binned*)malloc(M * binned_sbsbze(fold));
      if(beta == 0.0){
        memset(YI, 0, M * binned_sbsbze(fold));
      }else if(beta == 1.0){
        for(i = 0; i < M; i++){
          binned_sbsconv(fold, Y[i * incY], YI + i * binned_sbnum(fold));
        }
      }else{
        for(i = 0; i < M; i++){
          binned_sbsconv(fold, Y[i * incY] * beta, YI + i * binned_sbnum(fold));
        }
      }
      binnedBLAS_sbsgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < M; i++){
        Y[i * incY] = binned_ssbconv(fold, YI + i * binned_sbnum(fold));
      }
      break;
    default:
      YI = (float_binned*)malloc(N * binned_sbsbze(fold));
      if(beta == 0.0){
        memset(YI, 0, N * binned_sbsbze(fold));
      }else if(beta == 1.0){
        for(i = 0; i < N; i++){
          binned_sbsconv(fold, Y[i * incY], YI + i * binned_sbnum(fold));
        }
      }else{
        for(i = 0; i < N; i++){
          binned_sbsconv(fold, Y[i * incY] * beta, YI + i * binned_sbnum(fold));
        }
      }
      binnedBLAS_sbsgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < N; i++){
        Y[i * incY] = binned_ssbconv(fold, YI + i * binned_sbnum(fold));
      }
      break;
  }
  free(YI);
}
