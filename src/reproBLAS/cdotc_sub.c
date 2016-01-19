#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible conjugated dot product of complex single precision vectors X and Y
 *
 * Return the sum of the pairwise products of X and conjugated Y.
 *
 * The reproducible dot product is computed with indexed types of default fold using #idxdBLAS_cicdotc()
 *
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y complex single precision vector
 * @param incY Y vector stride (use every incY'th element)
 * @param dotc scalar return
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void reproBLAS_cdotc_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  reproBLAS_rcdotc_sub(SIDEFAULTFOLD, N, X, incX, Y, incY, dotc);
}
