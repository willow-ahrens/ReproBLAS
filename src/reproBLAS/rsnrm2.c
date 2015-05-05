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
  float_indexed *nrmi = sialloc(DEFAULT_FOLD);
  float scale;
  float nrm2;

  sisetzero(DEFAULT_FOLD, nrmi);

  scale = sisnrm(DEFAULT_FOLD, N, X, incX, nrmi);

  nrm2 = scale * sqrt(ssiconv(DEFAULT_FOLD, nrmi));
  free(nrmi);
  return nrm2;
}
