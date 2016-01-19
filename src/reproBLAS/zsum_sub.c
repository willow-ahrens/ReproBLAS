#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible sum of complex double precision vector X
 *
 * Return the sum of X.
 *
 * The reproducible sum is computed with indexed types of default fold using #idxdBLAS_zizsum()
 *
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param sum scalar return
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void reproBLAS_zsum_sub(const int N, const void* X, const int incX, void *sum) {
  reproBLAS_rzsum_sub(DIDEFAULTFOLD, N, X, incX, sum);
}
