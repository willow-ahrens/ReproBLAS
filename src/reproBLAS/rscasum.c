#include <reproBLAS.h>
#include <idxdBLAS.h>

#include "../../config.h"

float reproBLAS_rscasum(const int N, const void* X, const int incX) {
  float_indexed *asumi = idxd_sialloc(SIDEFAULTFOLD);
  float asum;

  idxd_sisetzero(SIDEFAULTFOLD, asumi);

  idxdBLAS_sicasum(SIDEFAULTFOLD, N, X, incX, asumi);

  asum = idxd_ssiconv(SIDEFAULTFOLD, asumi);
  free(asumi);
  return asum;
}

