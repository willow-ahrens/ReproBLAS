#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

float rscasum(const int N, const void* X, const int incX) {
  float_indexed *asumi = idxd_sialloc(SIDEFAULTFOLD);
  float asum;

  idxd_sisetzero(SIDEFAULTFOLD, asumi);

  sicasum(SIDEFAULTFOLD, N, X, incX, asumi);

  asum = idxd_ssiconv(SIDEFAULTFOLD, asumi);
  free(asumi);
  return asum;
}

