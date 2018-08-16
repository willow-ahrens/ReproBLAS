#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible unconjugated dot product of complex single precision vectors X and Y
 *
 * Return the sum of the pairwise products of X and Y.
 *
 * The reproducible dot product is computed with binned types using #binnedBLAS_cbcdotu()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y complex single precision vector
 * @param incY Y vector stride (use every incY'th element)
 * @param dotu scalar return
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void reproBLAS_rcdotu_sub(const int fold, const int N, const void* X, const int incX, const void *Y, const int incY, void *dotu) {
  float_complex_binned *dotui = binned_cballoc(fold);

  binned_cbsetzero(fold, dotui);

  binnedBLAS_cbcdotu(fold, N, X, incX, Y, incY, dotui);

  binned_ccbconv_sub(fold, dotui, dotu);
  free(dotui);
  return;
}
