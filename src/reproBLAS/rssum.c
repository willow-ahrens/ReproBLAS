#include <reproBLAS.h>
#include <idxdBLAS.h>

/**
 * @brief Compute the reproducible sum of single precision vector X
 *
 * Return the sum of X.
 *
 * The reproducible sum is computed with indexed types using #idxdBLAS_sissum()
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return sum of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
float reproBLAS_rssum(const int fold, const int N, const float* X, const int incX) {
  float_indexed *sumi = idxd_sialloc(fold);
  float sum;

  idxd_sisetzero(fold, sumi);

  idxdBLAS_sissum(fold, N, X, incX, sumi);

  sum = idxd_ssiconv(fold, sumi);
  free(sumi);
  return sum;
}
