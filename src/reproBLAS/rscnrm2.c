#include <math.h>

#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible Euclidian norm of complex single precision vector X
 *
 * Return the square root of the sum of the squared elements of X.
 *
 * The reproducible Euclidian norm is computed with scaled binned types using #binnedBLAS_sbcssq()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return Euclidian norm of X
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
float reproBLAS_rscnrm2(const int fold, const int N, const void* X, const int incX) {
  float_binned *ssq = binned_sballoc(fold);
  float scl;
  float nrm2;

  binned_sbsetzero(fold, ssq);

  scl = binnedBLAS_sbcssq(fold, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(binned_ssbconv(fold, ssq));
  free(ssq);
  return nrm2;
}
