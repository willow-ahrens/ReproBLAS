#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible sum of single precision vector X
 *
 * Return the sum of X.
 *
 * The reproducible sum is computed with binned types of default fold using #binnedBLAS_sbssum()
 *
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return sum of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
float reproBLAS_ssum(const int N, const float* X, const int incX) {
  return reproBLAS_rssum(SIDEFAULTFOLD, N, X, incX);
}
