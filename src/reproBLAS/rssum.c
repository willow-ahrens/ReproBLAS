/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

float rssum(const int N, const float* X, const int incX) {
  float_indexed *sumi = sialloc(DEFAULT_FOLD);
  float sum;

  sisetzero(DEFAULT_FOLD, sumi);

  sissum(DEFAULT_FOLD, N, X, incX, sumi);

  sum = ssiconv(DEFAULT_FOLD, sumi);
  free(sumi);
  return sum;
}

