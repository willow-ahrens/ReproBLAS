#include "binnedBLAS.h"

/**
 * @brief Add to binned double precision Y the absolute sum of double precision vector X
 *
 * Add to Y the binned sum of absolute values of elements in X.
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y binned scalar Y
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
void binnedBLAS_dbdasum(const int fold, const int N, const double *X, const int incX, double_binned *Y){
  binnedBLAS_dmdasum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
