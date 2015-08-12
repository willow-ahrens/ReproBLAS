#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

double rdzasum(const int N, const void* X, const int incX) {
  double_indexed *asumi = idxd_dialloc(DIDEFAULTFOLD);
  double asum;

  idxd_disetzero(DIDEFAULTFOLD, asumi);

  dizasum(DIDEFAULTFOLD, N, X, incX, asumi);

  asum = idxd_ddiconv(DIDEFAULTFOLD, asumi);
  free(asumi);
  return asum;
}

