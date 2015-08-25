/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "idxdBLAS.h"

void idxdBLAS_cicdotu(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float_complex_indexed *Z){
  idxdBLAS_cmcdotu(fold, N, X, incX, Y, incY, Z, 1, Z + 2 * fold, 1);
}
