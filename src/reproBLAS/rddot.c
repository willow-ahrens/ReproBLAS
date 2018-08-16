#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible dot product of double precision vectors X and Y
 *
 * Return the sum of the pairwise products of X and Y.
 *
 * The reproducible dot product is computed with binned types using #binnedBLAS_dbddot()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y double precision vector
 * @param incY Y vector stride (use every incY'th element)
 * @return the dot product of X and Y
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_rddot(const int fold, const int N, const double* X, const int incX, const double *Y, const int incY) {
  double_binned *doti = binned_dballoc(fold);
  double dot;

  binned_dbsetzero(fold, doti);

  binnedBLAS_dbddot(fold, N, X, incX, Y, incY, doti);

  dot = binned_ddbconv(fold, doti);
  free(doti);
  return dot;
}
