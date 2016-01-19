#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Add to indexed double precision matrix C the matrix-matrix product of double precision matrices A and B
 *
 * Performs one of the matrix-matrix operations
 *
 *   C := alpha*op(A)*op(B) + C,
 *
 * where  op(X) is one of
 *
 *   op(X) = X   or   op(X) = X**T,
 *
 * alpha is a scalar, A and B are matrices with op(A) an m by k matrix and op(B) a k by n matrix, and C is an indexed m by n matrix.
 *
 * @param fold the fold of the indexed types
 * @param Order a character specifying the matrix ordering ('r' or 'R' for row-major, 'c' or 'C' for column major)
 * @param TransA a character specifying whether or not to transpose A before taking the matrix-matrix product ('n' or 'N' not to transpose, 't' or 'T' or 'c' or 'C' to transpose)
 * @param TransB a character specifying whether or not to transpose B before taking the matrix-matrix product ('n' or 'N' not to transpose, 't' or 'T' or 'c' or 'C' to transpose)
 * @param M number of rows of matrix op(A) and of the matrix C.
 * @param N number of columns of matrix op(B) and of the matrix C.
 * @param K number of columns of matrix op(A) and columns of the matrix op(B).
 * @param alpha scalar alpha
 * @param A double precision matrix of dimension (ma, lda) in row-major or (lda, na) in column-major. (ma, na) is (M, K) if A is not transposed and (K, M) otherwise.
 * @param lda the first dimension of A as declared in the calling program. lda must be at least na in row major or ma in column major.
 * @param B double precision matrix of dimension (mb, ldb) in row-major or (ldb, nb) in column-major. (mb, nb) is (K, N) if B is not transposed and (N, K) otherwise.
 * @param ldb the first dimension of B as declared in the calling program. ldb must be at least nb in row major or mb in column major.
 * @param C indexed double precision matrix of dimension (M, ldc) in row-major or (ldc, N) in column-major.
 * @param ldc the first dimension of C as declared in the calling program. ldc must be at least N in row major or M in column major.
 *
 * @author Peter Ahrens
 * @date   18 Jan 2016
 */
void reproBLAS_dgemm(const char Order, const char TransA, const char TransB,
                     const int M, const int N, const int K,
                     const double alpha, const double *A, const int lda,
                     const double *B, const int ldb,
                     const double beta, double *C, const int ldc){
  reproBLAS_rdgemm(DIDEFAULTFOLD, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
