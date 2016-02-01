#include <reproBLAS.h>
#include <idxdBLAS.h>

/**
 * @brief Compute the reproducible absolute sum of single precision vector X
 *
 * Return the sum of absolute values of elements in X.
 *
 * The reproducible absolute sum is computed with indexed types using #idxdBLAS_sisasum()
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return absolute sum of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
float reproBLAS_rsasum(const int fold, const int N, const float* X, const int incX) {
  float_indexed *asumi = idxd_sialloc(fold);
  float asum;

  idxd_sisetzero(fold, asumi);

  idxdBLAS_sisasum(fold, N, X, incX, asumi);

  asum = idxd_ssiconv(fold, asumi);
  free(asumi);
  return asum;
}
