#include <reproBLAS.h>

#include "../../config.h"

/**
 * @brief Compute the reproducible unconjugated dot product of complex double precision vectors X and Y
 *
 * Return the sum of the pairwise products of X and Y.
 *
 * The reproducible dot product is computed with binned types of default fold using #binnedBLAS_zbzdotu()
 *
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y complex double precision vector
 * @param incY Y vector stride (use every incY'th element)
 * @param dotu scalar return
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
void reproBLAS_zdotu_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotu) {
  reproBLAS_rzdotu_sub(DIDEFAULTFOLD, N, X, incX, Y, incY, dotu);
}
