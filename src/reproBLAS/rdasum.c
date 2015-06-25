#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

double rdasum(const int N, const double* X, const int incX) {
  double_indexed *asumi = dialloc(DIDEFAULTFOLD);
  double asum;

  disetzero(DIDEFAULTFOLD, asumi);

  didasum(DIDEFAULTFOLD, N, X, incX, asumi);

  asum = ddiconv(DIDEFAULTFOLD, asumi);
  free(asumi);
  return asum;
}

