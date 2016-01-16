/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "idxdBLAS.h"

/**
 * @brief Compute indexed single precision sum Z of single precision vector X
 *
 * Set Z to the indexed sum of X.
 *
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Z indexed scalar Z
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void idxdBLAS_sissum(const int fold, const int N, const float *X, const int incX, float_indexed *Y){
  idxdBLAS_smssum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
