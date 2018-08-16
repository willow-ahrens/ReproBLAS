#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible absolute sum of double precision vector X
 *
 * Return the sum of absolute values of elements in X.
 *
 * The reproducible absolute sum is computed with binned types of default fold using #binnedBLAS_dbdasum()
 *
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return absolute sum of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_dasum(const int N, const double* X, const int incX) {
  return reproBLAS_rdasum(DIDEFAULTFOLD, N, X, incX);
}
