/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "idxdBLAS.h"

void idxdBLAS_sisdot(const int fold, const int N, const float *X, const int incX, const float *Y, const int incY, float_indexed *Z){
  idxdBLAS_smsdot(fold, N, X, incX, Y, incY, Z, 1, Z + fold, 1);
}
