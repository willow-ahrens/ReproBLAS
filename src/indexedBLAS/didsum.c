/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "indexedBLAS.h"

void didsum(const int fold, const int N, const double *X, const int incX, double_indexed *Y){
  dmdsum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
