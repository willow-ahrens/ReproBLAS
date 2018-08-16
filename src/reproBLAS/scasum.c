#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible absolute sum of complex single precision vector X
 *
 * Return the sum of magnitudes of elements of X.
 *
 * The reproducible absolute sum is computed with binned types of default fold using #binnedBLAS_sbcasum()
 *
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return absolute sum of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
float reproBLAS_scasum(const int N, const void* X, const int incX) {
  return reproBLAS_rscasum(SIDEFAULTFOLD, N, X, incX);
}
