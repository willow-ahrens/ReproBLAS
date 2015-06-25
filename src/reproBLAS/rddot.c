#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

double rddot(const int N, const double* X, const int incX, const double *Y, const int incY) {
  double_indexed *doti = dialloc(DIDEFAULTFOLD);
  double dot;

  disetzero(DIDEFAULTFOLD, doti);

  diddot(DIDEFAULTFOLD, N, X, incX, Y, incY, doti);

  dot = ddiconv(DIDEFAULTFOLD, doti);
  free(doti);
  return dot;
}

