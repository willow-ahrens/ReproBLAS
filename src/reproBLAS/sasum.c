#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible absolute sum of single precision vector X
 *
 * Return the sum of absolute values of elements in X.
 *
 * The reproducible absolute sum is computed with binned types of default fold using #binnedBLAS_sbsasum()
 *
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return absolute sum of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
float reproBLAS_sasum(const int N, const float* X, const int incX) {
  return reproBLAS_rsasum(SIDEFAULTFOLD, N, X, incX);
}
