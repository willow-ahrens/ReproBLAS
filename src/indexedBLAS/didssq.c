/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "indexedBLAS.h"

double didssq(const int fold, const int N, const double *X, const int incX, const double scaleY, double_indexed *Y){
  return dmdssq(fold, N, X, incX, scaleY, Y, 1, Y + fold, 1);
}
