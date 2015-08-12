/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "indexedBLAS.h"

void idxdBLAS_didasum(const int fold, const int N, const double *X, const int incX, double_indexed *Y){
  idxdBLAS_dmdasum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
