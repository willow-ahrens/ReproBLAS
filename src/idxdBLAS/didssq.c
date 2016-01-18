#include "idxdBLAS.h"

/**
 * @brief Add to scaled indexed double precision Y the scaled sum of squares of elements of double precision vector X
 *
 * Add to Y the scaled indexed sum of the squares of each element of X. The scaling of each square is performed using #idxd_dscale()
 *
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param scaleY the scaling factor of Y
 * @param Y indexed scalar Y
 * @return the new scaling factor of Y
 *
 * @author Peter Ahrens
 * @date   18 Jan 2016
 */
double idxdBLAS_didssq(const int fold, const int N, const double *X, const int incX, const double scaleY, double_indexed *Y){
  return idxdBLAS_dmdssq(fold, N, X, incX, scaleY, Y, 1, Y + fold, 1);
}
