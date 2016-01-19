#include <math.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

/**
 * @brief Compute the reproducible Euclidian norm of complex double precision vector X
 *
 * Return the square root of the sum of the squared elements of X.
 *
 * The reproducible Euclidian norm is computed with scaled indexed types using #idxdBLAS_dizssq()
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return Euclidian norm of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_rdznrm2(const int fold, const int N, const void* X, const int incX) {
  double_indexed *ssq = idxd_dialloc(fold);
  double scl;
  double nrm2;

  idxd_disetzero(fold, ssq);

  scl = idxdBLAS_dizssq(fold, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(idxd_ddiconv(fold, ssq));
  free(ssq);
  return nrm2;
}
