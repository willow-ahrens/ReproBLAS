#include <reproBLAS.h>
#include <idxdBLAS.h>

/**
 * @brief Compute the reproducible sum of complex double precision vector X
 *
 * Return the sum of X.
 *
 * The reproducible sum is computed with indexed types using #idxdBLAS_zizsum()
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param sum scalar return
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void reproBLAS_rzsum_sub(const int fold, const int N, const void* X, const int incX, void *sum) {
  double_complex_indexed *sumi = idxd_zialloc(fold);

  idxd_zisetzero(fold, sumi);

  idxdBLAS_zizsum(fold, N, X, incX, sumi);

  idxd_zziconv_sub(fold, sumi, sum);
  free(sumi);
  return;
}
