#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

double rdsum(const int N, const double* X, const int incX) {
  double_indexed *sumi = idxd_dialloc(DIDEFAULTFOLD);
  double sum;

  idxd_disetzero(DIDEFAULTFOLD, sumi);

  idxdBLAS_didsum(DIDEFAULTFOLD, N, X, incX, sumi);

  sum = idxd_ddiconv(DIDEFAULTFOLD, sumi);
  free(sumi);
  return sum;
}

