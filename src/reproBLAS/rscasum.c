#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible absolute sum of complex single precision vector X
 *
 * Return the sum of magnitudes of elements of X.
 *
 * The reproducible absolute sum is computed with binned types using #binnedBLAS_sbcasum()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return absolute sum of X
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
float reproBLAS_rscasum(const int fold, const int N, const void* X, const int incX) {
  float_binned *asumi = binned_sballoc(fold);
  float asum;

  binned_sbsetzero(fold, asumi);

  binnedBLAS_sbcasum(fold, N, X, incX, asumi);

  asum = binned_ssbconv(fold, asumi);
  free(asumi);
  return asum;
}
