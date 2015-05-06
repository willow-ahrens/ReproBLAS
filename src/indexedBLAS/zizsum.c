/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "indexedBLAS.h"

void zizsum(const int fold, const int N, const void *X, const int incX, double_complex_indexed *Y){
  zmzsum(fold, N, X, incX, Y, 1, Y + 2 * fold, 1);
}
