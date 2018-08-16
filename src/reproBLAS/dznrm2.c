#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible Euclidian norm of complex double precision vector X
 *
 * Return the square root of the sum of the squared elements of X.
 *
 * The reproducible Euclidian norm is computed with scaled binned types of default fold using #binnedBLAS_dbzssq()
 *
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return Euclidian norm of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_dznrm2(const int N, const void* X, const int incX) {
  return reproBLAS_rdznrm2(DIDEFAULTFOLD, N, X, incX);
}
