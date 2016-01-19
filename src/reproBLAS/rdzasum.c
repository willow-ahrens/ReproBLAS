#include <reproBLAS.h>
#include <idxdBLAS.h>

/**
 * @brief Compute the reproducible absolute sum of complex double precision vector X
 *
 * Return the sum of magnitudes of elements of X.
 *
 * The reproducible absolute sum is computed with indexed types using #idxdBLAS_dizasum()
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return absolute sum of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_rdzasum(const int fold, const int N, const void* X, const int incX) {
  double_indexed *asumi = idxd_dialloc(fold);
  double asum;

  idxd_disetzero(fold, asumi);

  idxdBLAS_dizasum(fold, N, X, incX, asumi);

  asum = idxd_ddiconv(fold, asumi);
  free(asumi);
  return asum;
}
