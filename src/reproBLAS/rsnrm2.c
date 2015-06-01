/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

float rsnrm2(const int N, const float* X, const int incX) {
  float_indexed *ssq = sialloc(DEFAULT_FOLD);
  float scl;
  float nrm2;

  sisetzero(DEFAULT_FOLD, ssq);

  scl = sisssq(DEFAULT_FOLD, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(ssiconv(DEFAULT_FOLD, ssq));
  free(ssq);
  return nrm2;
}
