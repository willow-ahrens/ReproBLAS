#include <reproBLAS.h>
#include <idxdBLAS.h>

/**
 * @brief Compute the reproducible sum of complex single precision vector X
 *
 * Return the sum of X.
 *
 * The reproducible sum is computed with indexed types using #idxdBLAS_cicsum()
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param sum scalar return
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void reproBLAS_rcsum_sub(const int fold, const int N, const void* X, const int incX, void *sum) {
  float_complex_indexed *sumi = idxd_cialloc(fold);

  idxd_cisetzero(fold, sumi);

  idxdBLAS_cicsum(fold, N, X, incX, sumi);

  idxd_cciconv_sub(fold, sumi, sum);
  free(sumi);
  return;
}
