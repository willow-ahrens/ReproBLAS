#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible sum of complex single precision vector X
 *
 * Return the sum of X.
 *
 * The reproducible sum is computed with binned types of default fold using #binnedBLAS_cbcsum()
 *
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param sum scalar return
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void reproBLAS_csum_sub(const int N, const void* X, const int incX, void *sum) {
  reproBLAS_rcsum_sub(SIDEFAULTFOLD, N, X, incX, sum);
}
