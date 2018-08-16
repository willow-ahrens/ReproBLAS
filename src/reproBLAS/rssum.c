#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible sum of single precision vector X
 *
 * Return the sum of X.
 *
 * The reproducible sum is computed with binned types using #binnedBLAS_sbssum()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return sum of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
float reproBLAS_rssum(const int fold, const int N, const float* X, const int incX) {
  float_binned *sumi = binned_sballoc(fold);
  float sum;

  binned_sbsetzero(fold, sumi);

  binnedBLAS_sbssum(fold, N, X, incX, sumi);

  sum = binned_ssbconv(fold, sumi);
  free(sumi);
  return sum;
}
