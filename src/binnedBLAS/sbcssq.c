#include "binnedBLAS.h"

/**
 * @brief Add to scaled binned single precision Y the scaled sum of squares of elements of complex single precision vector X
 *
 * Add to Y the scaled binned sum of the squares of each element of X. The scaling of each square is performed using #binned_sscale()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param scaleY the scaling factor of Y
 * @param Y binned scalar Y
 * @return the new scaling factor of Y
 *
 * @author Willow Ahrens
 * @date   18 Jan 2016
 */
float binnedBLAS_sbcssq(const int fold, const int N, const void *X, const int incX, const float scaleY, float_binned *Y){
  return binnedBLAS_smcssq(fold, N, X, incX, scaleY, Y, 1, Y + fold, 1);
}
