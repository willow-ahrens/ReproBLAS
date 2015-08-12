#include <stdlib.h>

#include <indexedBLAS.h>
#include <idxd.h>


void idxdBLAS_zizgemv(const int fold, const char Order,
             const char TransA, const int M, const int N,
             const void *alpha, const void *A, const int lda,
             const void *X, const int incX,
             double_complex_indexed *Y, const int incY){
  idxdBLAS_zmzgemv(fold, Order,
          TransA, M, N,
          alpha, A, lda,
          X, incX,
          Y, 1, 2 * incY,
          Y + 2 * fold, 1, 2 * incY);
}
