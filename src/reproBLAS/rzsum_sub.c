/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

void rzsum_sub(const int N, const void* X, const int incX, void *sum) {
  double_complex_indexed *sumi = zialloc(DEFAULT_FOLD);

  zisetzero(DEFAULT_FOLD, sumi);

  zizsum(DEFAULT_FOLD, N, X, incX, sumi);

  zziconv_sub(DEFAULT_FOLD, sumi, sum);
  free(sumi);
  return;
}

