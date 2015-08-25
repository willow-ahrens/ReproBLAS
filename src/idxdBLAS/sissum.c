/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "idxdBLAS.h"

void idxdBLAS_sissum(const int fold, const int N, const float *X, const int incX, float_indexed *Y){
  idxdBLAS_smssum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
