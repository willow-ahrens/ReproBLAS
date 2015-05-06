/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

void rcsum_sub(const int N, const void* X, const int incX, void *sum) {
  float_complex_indexed *sumi = cialloc(DEFAULT_FOLD);

  cisetzero(DEFAULT_FOLD, sumi);

  cicsum(DEFAULT_FOLD, N, X, incX, sumi);

  cciconv_sub(DEFAULT_FOLD, sumi, sum);
  free(sumi);
  return;
}

