#include "idxdBLAS.h"

/**
 * @brief Add to scaled indexed single precision Y the scaled sum of squares of elements of complex single precision vector X
 *
 * Add to Y the scaled indexed sum of the squares of each element of X. The scaling of each square is performed using #idxd_sscale()
 *
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param scaleY the scaling factor of Y
 * @param Y indexed scalar Y
 * @return the new scaling factor of Y
 *
 * @author Peter Ahrens
 * @date   18 Jan 2016
 */
float idxdBLAS_sicssq(const int fold, const int N, const void *X, const int incX, const float scaleY, float_indexed *Y){
  return idxdBLAS_smcssq(fold, N, X, incX, scaleY, Y, 1, Y + fold, 1);
}
