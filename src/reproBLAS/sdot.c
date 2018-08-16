#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible dot product of single precision vectors X and Y
 *
 * Return the sum of the pairwise products of X and Y.
 *
 * The reproducible dot product is computed with binned types of default fold using #binnedBLAS_sbsdot()
 *
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y single precision vector
 * @param incY Y vector stride (use every incY'th element)
 * @return the dot product of X and Y
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
float reproBLAS_sdot(const int N, const float* X, const int incX, const float *Y, const int incY) {
  return reproBLAS_rsdot(SIDEFAULTFOLD, N, X, incX, Y, incY);
}
