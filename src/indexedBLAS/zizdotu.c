/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "indexedBLAS.h"

void zizdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, double_complex_indexed *Z){
  zmzdotu(fold, N, X, incX, Y, incY, Z, 1, Z + 2 * fold, 1);
}
