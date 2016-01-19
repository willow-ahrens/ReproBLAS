#include "idxdBLAS.h"

/**
 * @brief Add to scaled indexed single precision Y the scaled sum of squares of elements of single precision vector X
 *
 * Add to Y the scaled indexed sum of the squares of each element of X. The scaling of each square is performed using #idxd_sscale()
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param scaleY the scaling factor of Y
 * @return the new scaling factor of Y
 *
 * @author Peter Ahrens
 * @date   18 Jan 2016
 */
float idxdBLAS_sisssq(const int fold, const int N, const float *X, const int incX, const float scaleY, float_indexed *Y){
  return idxdBLAS_smsssq(fold, N, X, incX, scaleY, Y, 1, Y + fold, 1);
}
