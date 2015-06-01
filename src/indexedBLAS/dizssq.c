/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "indexedBLAS.h"

double dizssq(const int fold, const int N, const void *X, const int incX, const double scaleY, double_indexed *Y){
  return dmzssq(fold, N, X, incX, scaleY, Y, 1, Y + fold, 1);
}
