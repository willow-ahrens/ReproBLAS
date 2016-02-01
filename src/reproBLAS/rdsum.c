#include <reproBLAS.h>
#include <idxdBLAS.h>

/**
 * @brief Compute the reproducible sum of double precision vector X
 *
 * Return the sum of X.
 *
 * The reproducible sum is computed with indexed types using #idxdBLAS_didsum()
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return sum of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_rdsum(const int fold, const int N, const double* X, const int incX) {
  double_indexed *sumi = idxd_dialloc(fold);
  double sum;

  idxd_disetzero(fold, sumi);

  idxdBLAS_didsum(fold, N, X, incX, sumi);

  sum = idxd_ddiconv(fold, sumi);
  free(sumi);
  return sum;
}
