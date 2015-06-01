/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "indexedBLAS.h"

float sisssq(const int fold, const int N, const float *X, const int incX, const float scaleY, float_indexed *Y){
  return smsssq(fold, N, X, incX, scaleY, Y, 1, Y + fold, 1);
}
