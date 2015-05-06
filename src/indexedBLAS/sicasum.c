/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "indexedBLAS.h"

void sicasum(const int fold, const int N, const void *X, const int incX, float_indexed *Y){
  smcasum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
