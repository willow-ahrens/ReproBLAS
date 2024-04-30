#include "binnedBLAS.h"

/**
 * @brief Add to binned complex double precision Z the conjugated dot product of complex double precision vectors X and Y
 *
 * Add to Z the binned sum of the pairwise products of X and conjugated Y.
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y complex double precision vector
 * @param incY Y vector stride (use every incY'th element)
 * @param Z scalar return Z
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
void binnedBLAS_zbzdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double_complex_binned *Z){
  binnedBLAS_zmzdotc(fold, N, X, incX, Y, incY, Z, 1, Z + 2 * fold, 1);
}
