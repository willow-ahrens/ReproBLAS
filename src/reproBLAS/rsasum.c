#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

float rsasum(const int N, const float* X, const int incX) {
  float_indexed *asumi = idxd_sialloc(SIDEFAULTFOLD);
  float asum;

  idxd_sisetzero(SIDEFAULTFOLD, asumi);

  sisasum(SIDEFAULTFOLD, N, X, incX, asumi);

  asum = idxd_ssiconv(SIDEFAULTFOLD, asumi);
  free(asumi);
  return asum;
}

