/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "idxdBLAS.h"

double idxdBLAS_dizssq(const int fold, const int N, const void *X, const int incX, const double scaleY, double_indexed *Y){
  return idxdBLAS_dmzssq(fold, N, X, incX, scaleY, Y, 1, Y + fold, 1);
}
