#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible sum of double precision vector X
 *
 * Return the sum of X.
 *
 * The reproducible sum is computed with binned types using #binnedBLAS_dbdsum()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return sum of X
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_rdsum(const int fold, const int N, const double* X, const int incX) {
  double_binned *sumi = binned_dballoc(fold);
  double sum;

  binned_dbsetzero(fold, sumi);

  binnedBLAS_dbdsum(fold, N, X, incX, sumi);

  sum = binned_ddbconv(fold, sumi);
  free(sumi);
  return sum;
}
