#include <reproBLAS.h>
#include <binnedBLAS.h>

/**
 * @brief Compute the reproducible sum of complex double precision vector X
 *
 * Return the sum of X.
 *
 * The reproducible sum is computed with binned types using #binnedBLAS_zbzsum()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param sum scalar return
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
void reproBLAS_rzsum_sub(const int fold, const int N, const void* X, const int incX, void *sum) {
  double_complex_binned *sumi = binned_zballoc(fold);

  binned_zbsetzero(fold, sumi);

  binnedBLAS_zbzsum(fold, N, X, incX, sumi);

  binned_zzbconv_sub(fold, sumi, sum);
  free(sumi);
  return;
}
