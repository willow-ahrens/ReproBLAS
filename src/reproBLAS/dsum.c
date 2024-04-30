#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible sum of double precision vector X
 *
 * Return the sum of X.
 *
 * The reproducible sum is computed with binned types of default fold using #binnedBLAS_dbdsum()
 *
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return sum of X
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_dsum(const int N, const double* X, const int incX) {
  return reproBLAS_rdsum(DIDEFAULTFOLD, N, X, incX);
}
