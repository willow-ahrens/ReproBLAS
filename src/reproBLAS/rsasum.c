#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible absolute sum of single precision vector X
 *
 * Return the sum of absolute values of elements in X.
 *
 * The reproducible absolute sum is computed with binned types using #binnedBLAS_sbsasum()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return absolute sum of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
float reproBLAS_rsasum(const int fold, const int N, const float* X, const int incX) {
  float_binned *asumi = binned_sballoc(fold);
  float asum;

  binned_sbsetzero(fold, asumi);

  binnedBLAS_sbsasum(fold, N, X, incX, asumi);

  asum = binned_ssbconv(fold, asumi);
  free(asumi);
  return asum;
}
