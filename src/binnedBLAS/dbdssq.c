#include "binnedBLAS.h"

/**
 * @brief Add to scaled binned double precision Y the scaled sum of squares of elements of double precision vector X
 *
 * Add to Y the scaled binned sum of the squares of each element of X. The scaling of each square is performed using #binned_dscale()
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param scaleY the scaling factor of Y
 * @param Y binned scalar Y
 * @return the new scaling factor of Y
 *
 * @author Willow Ahrens
 * @date   18 Jan 2016
 */
double binnedBLAS_dbdssq(const int fold, const int N, const double *X, const int incX, const double scaleY, double_binned *Y){
  return binnedBLAS_dmdssq(fold, N, X, incX, scaleY, Y, 1, Y + fold, 1);
}
