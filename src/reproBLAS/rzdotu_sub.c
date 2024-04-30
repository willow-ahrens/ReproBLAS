#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible unconjugated dot product of complex double precision vectors X and Y
 *
 * Return the sum of the pairwise products of X and Y.
 *
 * The reproducible dot product is computed with binned types using #binnedBLAS_zbzdotu()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y complex double precision vector
 * @param incY Y vector stride (use every incY'th element)
 * @param dotu scalar return
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
void reproBLAS_rzdotu_sub(const int fold, const int N, const void* X, const int incX, const void *Y, const int incY, void *dotu) {
  double_complex_binned *dotui = binned_zballoc(fold);

  binned_zbsetzero(fold, dotui);

  binnedBLAS_zbzdotu(fold, N, X, incX, Y, incY, dotui);

  binned_zzbconv_sub(fold, dotui, dotu);
  free(dotui);
  return;
}
