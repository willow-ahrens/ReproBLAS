/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "idxdBLAS.h"

double idxdBLAS_didssq(const int fold, const int N, const double *X, const int incX, const double scaleY, double_indexed *Y){
  return idxdBLAS_dmdssq(fold, N, X, incX, scaleY, Y, 1, Y + fold, 1);
}
