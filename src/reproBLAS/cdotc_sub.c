#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute reproducible conjugated dot product of complex single precision vectors X and Y
 *
 * Returns the reproducible sum of the pairwise products of X and conjugated Y.
 *
 * The dot product is computed using indexed types of default fold with #idxdBLAS_cicdotu().
 *
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y complex single precision vector
 * @param incY Y vector stride (use every incY'th element)
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void reproBLAS_cdotc_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  reproBLAS_rcdotc_sub(SIDEFAULTFOLD, N, X, incX, Y, incY, dotc);
}
