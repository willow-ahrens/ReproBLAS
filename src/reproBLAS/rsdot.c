#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible dot product of single precision vectors X and Y
 *
 * Return the sum of the pairwise products of X and Y.
 *
 * The reproducible dot product is computed with binned types using #binnedBLAS_sbsdot()
 *
 * @param fold the fold of the binned types
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
float reproBLAS_rsdot(const int fold, const int N, const float* X, const int incX, const float *Y, const int incY) {
  float_binned *doti = binned_sballoc(fold);
  float dot;

  binned_sbsetzero(fold, doti);

  binnedBLAS_sbsdot(fold, N, X, incX, Y, incY, doti);

  dot = binned_ssbconv(fold, doti);
  free(doti);
  return dot;
}
