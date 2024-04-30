#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible absolute sum of complex double precision vector X
 *
 * Return the sum of magnitudes of elements of X.
 *
 * The reproducible absolute sum is computed with binned types using #binnedBLAS_dbzasum()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return absolute sum of X
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_rdzasum(const int fold, const int N, const void* X, const int incX) {
  double_binned *asumi = binned_dballoc(fold);
  double asum;

  binned_dbsetzero(fold, asumi);

  binnedBLAS_dbzasum(fold, N, X, incX, asumi);

  asum = binned_ddbconv(fold, asumi);
  free(asumi);
  return asum;
}
