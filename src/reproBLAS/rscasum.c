#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

float rscasum(const int N, const void* X, const int incX) {
  float_indexed *asumi = sialloc(SIDEFAULTFOLD);
  float asum;

  sisetzero(SIDEFAULTFOLD, asumi);

  sicasum(SIDEFAULTFOLD, N, X, incX, asumi);

  asum = ssiconv(SIDEFAULTFOLD, asumi);
  free(asumi);
  return asum;
}

