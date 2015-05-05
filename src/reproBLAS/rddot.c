/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

double rddot(const int N, const double* X, const int incX, const double *Y, const int incY) {
  double_indexed *doti = dialloc(DEFAULT_FOLD);
  double dot;

  disetzero(DEFAULT_FOLD, doti);

  diddot(DEFAULT_FOLD, N, X, incX, Y, incY, doti);

  dot = ddiconv(DEFAULT_FOLD, doti);
  free(doti);
  return dot;
}

