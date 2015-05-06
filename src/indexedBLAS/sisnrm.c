/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "indexedBLAS.h"

float sisnrm(const int fold, const int N, const float *X, const int incX, float_indexed *Y){
  return smsnrm(fold, N, X, incX, Y, 1, Y + fold, 1);
}
