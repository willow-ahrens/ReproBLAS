#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible conjugated dot product of complex single precision vectors X and Y
 *
 * Return the sum of the pairwise products of X and conjugated Y.
 *
 * The reproducible dot product is computed with binned types using #binnedBLAS_cbcdotc()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y complex single precision vector
 * @param incY Y vector stride (use every incY'th element)
 * @param dotc scalar return
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
void reproBLAS_rcdotc_sub(const int fold, const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  float_complex_binned *dotci = binned_cballoc(fold);

  binned_cbsetzero(fold, dotci);

  binnedBLAS_cbcdotc(fold, N, X, incX, Y, incY, dotci);

  binned_ccbconv_sub(fold, dotci, dotc);
  free(dotci);
  return;
}
