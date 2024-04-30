#include "binnedBLAS.h"

/**
 * @brief Add to binned complex double precision Y the sum of complex double precision vector X
 *
 * Add to Y the binned sum of X.
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y binned scalar Y
 *
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
void binnedBLAS_zbzsum(const int fold, const int N, const void *X, const int incX, double_complex_binned *Y){
  binnedBLAS_zmzsum(fold, N, X, incX, Y, 1, Y + 2 * fold, 1);
}
