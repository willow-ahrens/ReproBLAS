/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

double rdnrm2(const int N, const double* X, const int incX) {
  double_indexed *ssq = dialloc(DEFAULT_FOLD);
  double scl;
  double nrm2;

  disetzero(DEFAULT_FOLD, ssq);

  scl = didssq(DEFAULT_FOLD, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(ddiconv(DEFAULT_FOLD, ssq));
  free(ssq);
  return nrm2;
}
