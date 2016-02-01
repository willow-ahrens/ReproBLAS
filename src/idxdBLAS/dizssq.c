#include "idxdBLAS.h"

/**
 * @brief Add to scaled indexed double precision Y the scaled sum of squares of elements of complex double precision vector X
 *
 * Add to Y the scaled indexed sum of the squares of each element of X. The scaling of each square is performed using #idxd_dscale()
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param scaleY the scaling factor of Y
 * @param Y indexed scalar Y
 * @return the new scaling factor of Y
 *
 * @author Peter Ahrens
 * @date   18 Jan 2016
 */
double idxdBLAS_dizssq(const int fold, const int N, const void *X, const int incX, const double scaleY, double_indexed *Y){
  return idxdBLAS_dmzssq(fold, N, X, incX, scaleY, Y, 1, Y + fold, 1);
}
