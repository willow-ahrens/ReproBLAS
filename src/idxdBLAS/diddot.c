/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "idxdBLAS.h"

void idxdBLAS_diddot(const int fold, const int N, const double *X, const int incX, const double *Y, const int incY, double_indexed *Z){
  idxdBLAS_dmddot(fold, N, X, incX, Y, incY, Z, 1, Z + fold, 1);
}
