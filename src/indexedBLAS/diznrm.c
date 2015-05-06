/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "indexedBLAS.h"

double diznrm(const int fold, const int N, const void *X, const int incX, double_indexed *Y){
  return dmznrm(fold, N, X, incX, Y, 1, Y + fold, 1);
}
