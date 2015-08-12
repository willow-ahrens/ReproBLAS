#include <math.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

#include "../../config.h"

double rdznrm2(const int N, const void* X, const int incX) {
  double_indexed *ssq = idxd_dialloc(DIDEFAULTFOLD);
  double scl;
  double nrm2;

  idxd_disetzero(DIDEFAULTFOLD, ssq);

  scl = idxdBLAS_dizssq(DIDEFAULTFOLD, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(idxd_ddiconv(DIDEFAULTFOLD, ssq));
  free(ssq);
  return nrm2;
}
