#include <math.h>

#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible Euclidian norm of complex double precision vector X
 *
 * Return the square root of the sum of the squared elements of X.
 *
 * The reproducible Euclidian norm is computed with scaled binned types using #binnedBLAS_dbzssq()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return Euclidian norm of X
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_rdznrm2(const int fold, const int N, const void* X, const int incX) {
  double_binned *ssq = binned_dballoc(fold);
  double scl;
  double nrm2;

  binned_dbsetzero(fold, ssq);

  scl = binnedBLAS_dbzssq(fold, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(binned_ddbconv(fold, ssq));
  free(ssq);
  return nrm2;
}
