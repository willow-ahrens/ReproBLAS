#include "binnedBLAS.h"

/**
 * @brief Add to binned double precision Z the dot product of double precision vectors X and Y
 *
 * Add to Z the binned sum of the pairwise products of X and Y.
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y double precision vector
 * @param incY Y vector stride (use every incY'th element)
 * @param Z binned scalar Z
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
void binnedBLAS_dbddot(const int fold, const int N, const double *X, const int incX, const double *Y, const int incY, double_binned *Z){
  binnedBLAS_dmddot(fold, N, X, incX, Y, incY, Z, 1, Z + fold, 1);
}
