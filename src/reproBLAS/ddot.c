#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible dot product of double precision vectors X and Y
 *
 * Return the sum of the pairwise products of X and Y.
 *
 * The reproducible dot product is computed with indexed types of default fold using #idxdBLAS_diddot()
 *
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y double precision vector
 * @param incY Y vector stride (use every incY'th element)
 * @return the dot product of X and Y
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
double reproBLAS_ddot(const int N, const double* X, const int incX, const double *Y, const int incY) {
  return reproBLAS_rddot(DIDEFAULTFOLD, N, X, incX, Y, incY);
}
