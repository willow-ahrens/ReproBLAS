/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

void rcdotc_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  float_complex_indexed *dotci = cialloc(DEFAULT_FOLD);

  cisetzero(DEFAULT_FOLD, dotci);

  cicdotc(DEFAULT_FOLD, N, X, incX, Y, incY, dotci);

  cciconv_sub(DEFAULT_FOLD, dotci, dotc);
  free(dotci);
  return;
}

