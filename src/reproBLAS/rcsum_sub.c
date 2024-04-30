#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible sum of complex single precision vector X
 *
 * Return the sum of X.
 *
 * The reproducible sum is computed with binned types using #binnedBLAS_cbcsum()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param sum scalar return
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
void reproBLAS_rcsum_sub(const int fold, const int N, const void* X, const int incX, void *sum) {
  float_complex_binned *sumi = binned_cballoc(fold);

  binned_cbsetzero(fold, sumi);

  binnedBLAS_cbcsum(fold, N, X, incX, sumi);

  binned_ccbconv_sub(fold, sumi, sum);
  free(sumi);
  return;
}
