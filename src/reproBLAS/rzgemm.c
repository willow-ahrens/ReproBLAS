#include <stdlib.h>
#include <string.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

/**
 * @brief Add to complex double precision matrix C the reproducible matrix-matrix product of complex double precision matrices A and B
 *
 * Performs one of the matrix-matrix operations
 *
 *   C := alpha*op(A)*op(B) + beta*C,
 *
 * where  op(X) is one of
 *
 *   op(X) = X   or   op(X) = X**T   or   op(X) = X**H,
 *
 * alpha and beta are scalars, A and B and C are matrices with op(A) an M by K matrix, op(B) a K by N matrix, and C is an M by N matrix.
 *
 * The matrix-matrix product is computed using indexed types with #idxdBLAS_zizgemm()
 *
 * @param fold the fold of the indexed types
 * @param Order a character specifying the matrix ordering ('r' or 'R' for row-major, 'c' or 'C' for column major)
 * @param TransA a character specifying whether or not to transpose A before taking the matrix-matrix product ('n' or 'N' not to transpose, 't' or 'T' to transpose, 'c' or 'C' to conjugate transpose)
 * @param TransB a character specifying whether or not to transpose B before taking the matrix-matrix product ('n' or 'N' not to transpose, 't' or 'T' to transpose, 'c' or 'C' to conjugate transpose)
 * @param M number of rows of matrix op(A) and of the matrix C.
 * @param N number of columns of matrix op(B) and of the matrix C.
 * @param K number of columns of matrix op(A) and columns of the matrix op(B).
 * @param alpha scalar alpha
 * @param A complex double precision matrix of dimension (ma, lda) in row-major or (lda, na) in column-major. (ma, na) is (M, K) if A is not transposed and (K, M) otherwise.
 * @param lda the first dimension of A as declared in the calling program. lda must be at least na in row major or ma in column major.
 * @param B complex double precision matrix of dimension (mb, ldb) in row-major or (ldb, nb) in column-major. (mb, nb) is (K, N) if B is not transposed and (N, K) otherwise.
 * @param ldb the first dimension of B as declared in the calling program. ldb must be at least nb in row major or mb in column major.
 * @param C complex double precision matrix of dimension (M, ldc) in row-major or (ldc, N) in column-major.
 * @param ldc the first dimension of C as declared in the calling program. ldc must be at least N in row major or M in column major.
 *
 * @author Peter Ahrens
 * @date   18 Jan 2016
 */
void reproBLAS_rzgemm(const int fold, const char Order, const char TransA, const char TransB,
                      const int M, const int N, const int K,
                      const void *alpha, const void *A, const int lda,
                      const void *B, const int ldb,
                      const void *beta, void *C, const int ldc){
  double_complex_indexed *CI;
  double betaC[2];
  int i;
  int j;

  if(M == 0 || N == 0){
    return;
  }

  CI = (double_complex_indexed*)malloc(M * N * idxd_zisize(fold));
  switch(Order){
    case 'r':
    case 'R':
      if(((double*)beta)[0] == 0.0 && ((double*)beta)[1] == 0.0){
        memset(CI, 0, M * N * idxd_zisize(fold));
      }else if(((double*)beta)[0] == 1.0 && ((double*)beta)[1] == 0.0){
        for(i = 0; i < M; i++){
          for(j = 0; j < N; j++){
            idxd_zizconv(fold, ((double*)C) + 2 * (i * ldc + j), CI + (i * N + j) * idxd_zinum(fold));
          }
        }
      }else{
        for(i = 0; i < M; i++){
          for(j = 0; j < N; j++){
            betaC[0] = ((double*)C)[2 * (i * ldc + j)] * ((double*)beta)[0] - ((double*)C)[2 * (i * ldc + j) + 1] * ((double*)beta)[1],
            betaC[1] = ((double*)C)[2 * (i * ldc + j)] * ((double*)beta)[1] + ((double*)C)[2 * (i * ldc + j) + 1] * ((double*)beta)[0],
            idxd_zizconv(fold, betaC, CI + (i * N + j) * idxd_zinum(fold));
          }
        }
      }
      idxdBLAS_zizgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, N);
      for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
          idxd_zziconv_sub(fold, CI + (i * N + j) * idxd_zinum(fold), ((double*)C) + 2 * (i * ldc + j));
        }
      }
      break;
    default:
      if(((double*)beta)[0] == 0.0 && ((double*)beta)[1] == 0.0){
        memset(CI, 0, M * N * idxd_zisize(fold));
      }else if(((double*)beta)[0] == 1.0 && ((double*)beta)[1] == 0.0){
        for(j = 0; j < N; j++){
          for(i = 0; i < M; i++){
            idxd_zizconv(fold, ((double*)C) + 2 * (j * ldc + i), CI + (j * M + i) * idxd_zinum(fold));
          }
        }
      }else{
        for(j = 0; j < N; j++){
          for(i = 0; i < M; i++){
            betaC[0] = ((double*)C)[2 * (j * ldc + i)] * ((double*)beta)[0] - ((double*)C)[2 * (j * ldc + i) + 1] * ((double*)beta)[1],
            betaC[1] = ((double*)C)[2 * (j * ldc + i)] * ((double*)beta)[1] + ((double*)C)[2 * (j * ldc + i) + 1] * ((double*)beta)[0],
            idxd_zizconv(fold, betaC, CI + (j * M + i) * idxd_zinum(fold));
          }
        }
      }
      idxdBLAS_zizgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, M);
      for(j = 0; j < N; j++){
        for(i = 0; i < M; i++){
          idxd_zziconv_sub(fold, CI + (j * M + i) * idxd_zinum(fold), ((double*)C) + 2 * (j * ldc + i));
        }
      }
      break;
  }
  free(CI);
}
