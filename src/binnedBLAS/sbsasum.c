#include "binnedBLAS.h"

/**
 * @brief Add to binned single precision Y the absolute sum of single precision vector X
 *
 * Add to Y the binned sum of absolute values of elements in X.
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y binned scalar Y
 * @author Willow Ahrens
 * @date   15 Jan 2016
 */
void binnedBLAS_sbsasum(const int fold, const int N, const float *X, const int incX, float_binned *Y){
  binnedBLAS_smsasum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
