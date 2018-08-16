#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible Euclidian norm of double precision vector X
 *
 * Return the square root of the sum of the squared elements of X.
 *
 * The reproducible Euclidian norm is computed with scaled binned types of default fold using #binnedBLAS_dbdssq()
 *
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return Euclidian norm of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_dnrm2(const int N, const double* X, const int incX) {
  return reproBLAS_rdnrm2(DIDEFAULTFOLD, N, X, incX);
}
