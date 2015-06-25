#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

double rdzasum(const int N, const void* X, const int incX) {
  double_indexed *asumi = dialloc(DIDEFAULTFOLD);
  double asum;

  disetzero(DIDEFAULTFOLD, asumi);

  dizasum(DIDEFAULTFOLD, N, X, incX, asumi);

  asum = ddiconv(DIDEFAULTFOLD, asumi);
  free(asumi);
  return asum;
}

