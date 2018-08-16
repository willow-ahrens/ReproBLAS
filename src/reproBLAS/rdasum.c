#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible absolute sum of double precision vector X
 *
 * Return the sum of absolute values of elements in X.
 *
 * The reproducible absolute sum is computed with binned types using #binnedBLAS_dbdasum()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return absolute sum of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_rdasum(const int fold, const int N, const double* X, const int incX) {
  double_binned *asumi = binned_dballoc(fold);
  double asum;

  binned_dbsetzero(fold, asumi);

  binnedBLAS_dbdasum(fold, N, X, incX, asumi);

  asum = binned_ddbconv(fold, asumi);
  free(asumi);
  return asum;
}
