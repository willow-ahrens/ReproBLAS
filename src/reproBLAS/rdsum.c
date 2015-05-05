/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

double rdsum(const int N, const double* X, const int incX) {
  double_indexed *sumi = dialloc(DEFAULT_FOLD);
  double sum;

  disetzero(DEFAULT_FOLD, sumi);

  didsum(DEFAULT_FOLD, N, X, incX, sumi);

  sum = ddiconv(DEFAULT_FOLD, sumi);
  free(sumi);
  return sum;
}

