#include <math.h>

#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

double rdznrm2(const int N, const void* X, const int incX) {
  double_indexed *ssq = dialloc(DIDEFAULTFOLD);
  double scl;
  double nrm2;

  disetzero(DIDEFAULTFOLD, ssq);

  scl = dizssq(DIDEFAULTFOLD, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(ddiconv(DIDEFAULTFOLD, ssq));
  free(ssq);
  return nrm2;
}
