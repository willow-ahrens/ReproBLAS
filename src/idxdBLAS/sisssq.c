/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "idxdBLAS.h"

float idxdBLAS_sisssq(const int fold, const int N, const float *X, const int incX, const float scaleY, float_indexed *Y){
  return idxdBLAS_smsssq(fold, N, X, incX, scaleY, Y, 1, Y + fold, 1);
}
