#include <reproBLAS.h>
#include <idxdBLAS.h>

#include "../../config.h"

double rddot(const int N, const double* X, const int incX, const double *Y, const int incY) {
  double_indexed *doti = idxd_dialloc(DIDEFAULTFOLD);
  double dot;

  idxd_disetzero(DIDEFAULTFOLD, doti);

  idxdBLAS_diddot(DIDEFAULTFOLD, N, X, incX, Y, incY, doti);

  dot = idxd_ddiconv(DIDEFAULTFOLD, doti);
  free(doti);
  return dot;
}

