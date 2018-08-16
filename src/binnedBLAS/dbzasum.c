#include "binnedBLAS.h"

/**
 * @brief Add to binned double precision Y the absolute sum of complex double precision vector X
 *
 * Add to Y the binned sum of magnitudes of elements of X.
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y binned scalar Y
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void binnedBLAS_dbzasum(const int fold, const int N, const void *X, const int incX, double_binned *Y){
  binnedBLAS_dmzasum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
