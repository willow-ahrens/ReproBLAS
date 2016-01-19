#include <reproBLAS.h>
#include <idxdBLAS.h>

/**
 * @brief Compute the reproducible conjugated dot product of complex single precision vectors X and Y
 *
 * Return the sum of the pairwise products of X and conjugated Y.
 *
 * The reproducible dot product is computed with indexed types using #idxdBLAS_cicdotc()
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y complex single precision vector
 * @param incY Y vector stride (use every incY'th element)
 * @param dotc scalar return
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void reproBLAS_rcdotc_sub(const int fold, const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  float_complex_indexed *dotci = idxd_cialloc(fold);

  idxd_cisetzero(fold, dotci);

  idxdBLAS_cicdotc(fold, N, X, incX, Y, incY, dotci);

  idxd_cciconv_sub(fold, dotci, dotc);
  free(dotci);
  return;
}
