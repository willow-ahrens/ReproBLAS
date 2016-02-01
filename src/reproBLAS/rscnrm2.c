#include <math.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

/**
 * @brief Compute the reproducible Euclidian norm of complex single precision vector X
 *
 * Return the square root of the sum of the squared elements of X.
 *
 * The reproducible Euclidian norm is computed with scaled indexed types using #idxdBLAS_sicssq()
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return Euclidian norm of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
float reproBLAS_rscnrm2(const int fold, const int N, const void* X, const int incX) {
  float_indexed *ssq = idxd_sialloc(fold);
  float scl;
  float nrm2;

  idxd_sisetzero(fold, ssq);

  scl = idxdBLAS_sicssq(fold, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(idxd_ssiconv(fold, ssq));
  free(ssq);
  return nrm2;
}
