#include <stdlib.h>
#include <string.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

/**
 * @brief Add to complex single precision vector Y the reproducible matrix-vector product of complex single precision matrix A and complex single precision vector X
 *
 * Performs one of the matrix-vector operations
 *
 *   y := alpha*A*x + beta*y   or   y := alpha*A**T*x + beta*y   or   y := alpha*A**H*x + beta*y,
 *
 * where alpha and beta are scalars, x and y are vectors, and A is an M by N matrix.
 *
 * The matrix-vector product is computed using indexed types with #idxdBLAS_cicgemv()
 *
 * @param fold the fold of the indexed types
 * @param Order a character specifying the matrix ordering ('r' or 'R' for row-major, 'c' or 'C' for column major)
 * @param TransA a character specifying whether or not to transpose A before taking the matrix-vector product ('n' or 'N' not to transpose, 't' or 'T' to transpose, 'c' or 'C' to conjugate transpose)
 * @param M number of rows of matrix A
 * @param N number of columns of matrix A
 * @param alpha scalar alpha
 * @param A complex single precision matrix of dimension (M, lda) in row-major or (lda, N) in column-major
 * @param lda the first dimension of A as declared in the calling program
 * @param X complex single precision vector of at least size N if not transposed or size M otherwise
 * @param incX X vector stride (use every incX'th element)
 * @param beta scalar beta
 * @param Y complex single precision vector Y of at least size M if not transposed or size N otherwise
 * @param incY Y vector stride (use every incY'th element)
 *
 * @author Peter Ahrens
 * @date   18 Jan 2016
 */
void reproBLAS_rcgemv(const int fold, const char Order,
                      const char TransA, const int M, const int N,
                      const void *alpha, const void *A, const int lda,
                      const void *X, const int incX,
                      const void *beta, void *Y, const int incY){
  float_complex_indexed *YI;
  float betaY[2];
  int i;

  if(N == 0 || M == 0){
    return;
  }

  switch(TransA){
    case 'n':
    case 'N':
      YI = (float_complex_indexed*)malloc(M * idxd_cisize(fold));
      if(((float*)beta)[0] == 0.0 && ((float*)beta)[1] == 0.0){
        memset(YI, 0, M * idxd_cisize(fold));
      }else if(((float*)beta)[0] == 1.0 && ((float*)beta)[1] == 0.0){
        for(i = 0; i < M; i++){
          idxd_cicconv(fold, ((float*)Y) + 2 * i * incY, YI + i * idxd_cinum(fold));
        }
      }else{
        for(i = 0; i < M; i++){
          betaY[0] = ((float*)Y)[2 * i * incY] * ((float*)beta)[0] - ((float*)Y)[2 * i * incY + 1] * ((float*)beta)[1];
          betaY[1] = ((float*)Y)[2 * i * incY] * ((float*)beta)[1] + ((float*)Y)[2 * i * incY + 1] * ((float*)beta)[0];
          idxd_cicconv(fold, betaY, YI + i * idxd_cinum(fold));
        }
      }
      idxdBLAS_cicgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < M; i++){
        idxd_cciconv_sub(fold, YI + i * idxd_cinum(fold), ((float*)Y) + 2 * i * incY);
      }
      break;
    default:
      YI = (float_complex_indexed*)malloc(N * idxd_cisize(fold));
      if(((float*)beta)[0] == 0.0 && ((float*)beta)[1] == 0.0){
        memset(YI, 0, N * idxd_cisize(fold));
      }else if(((float*)beta)[0] == 1.0 && ((float*)beta)[1] == 0.0){
        for(i = 0; i < N; i++){
          idxd_cicconv(fold, ((float*)Y) + 2 * i * incY, YI + i * idxd_cinum(fold));
        }
      }else{
        for(i = 0; i < N; i++){
          betaY[0] = ((float*)Y)[2 * i * incY] * ((float*)beta)[0] - ((float*)Y)[2 * i * incY + 1] * ((float*)beta)[1];
          betaY[1] = ((float*)Y)[2 * i * incY] * ((float*)beta)[1] + ((float*)Y)[2 * i * incY + 1] * ((float*)beta)[0];
          idxd_cicconv(fold, betaY, YI + i * idxd_cinum(fold));
        }
      }
      idxdBLAS_cicgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(i = 0; i < N; i++){
        idxd_cciconv_sub(fold, YI + i * idxd_cinum(fold), ((float*)Y) + 2 * i * incY);
      }
      break;
  }

  free(YI);
}
