/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

double rdznrm2(const int N, const void* X, const int incX) {
  double_indexed *nrmi = dialloc(DEFAULT_FOLD);
  double scale;
  double nrm2;

  disetzero(DEFAULT_FOLD, nrmi);

  scale = diznrm(DEFAULT_FOLD, N, X, incX, nrmi);

  nrm2 = scale * sqrt(ddiconv(DEFAULT_FOLD, nrmi));
  free(nrmi);
  return nrm2;
}
