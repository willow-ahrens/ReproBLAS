/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "indexedBLAS.h"

float idxdBLAS_sicssq(const int fold, const int N, const void *X, const int incX, const float scaleY, float_indexed *Y){
  return idxdBLAS_smcssq(fold, N, X, incX, scaleY, Y, 1, Y + fold, 1);
}
