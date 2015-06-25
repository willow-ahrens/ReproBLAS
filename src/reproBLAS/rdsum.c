#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

double rdsum(const int N, const double* X, const int incX) {
  double_indexed *sumi = dialloc(DIDEFAULTFOLD);
  double sum;

  disetzero(DIDEFAULTFOLD, sumi);

  didsum(DIDEFAULTFOLD, N, X, incX, sumi);

  sum = ddiconv(DIDEFAULTFOLD, sumi);
  free(sumi);
  return sum;
}

