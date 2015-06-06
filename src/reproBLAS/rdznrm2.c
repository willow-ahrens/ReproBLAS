#include <math.h>

#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

double rdznrm2(const int N, const void* X, const int incX) {
  double_indexed *ssq = dialloc(DEFAULT_FOLD);
  double scl;
  double nrm2;

  disetzero(DEFAULT_FOLD, ssq);

  scl = dizssq(DEFAULT_FOLD, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(ddiconv(DEFAULT_FOLD, ssq));
  free(ssq);
  return nrm2;
}
