#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible absolute sum of complex double precision vector X
 *
 * Return the sum of magnitudes of elements of X.
 *
 * The reproducible absolute sum is computed with indexed types of default fold using #idxdBLAS_dizasum()
 *
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @return absolute sum of X
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_dzasum(const int N, const void* X, const int incX) {
  return reproBLAS_rdzasum(DIDEFAULTFOLD, N, X, incX);
}
