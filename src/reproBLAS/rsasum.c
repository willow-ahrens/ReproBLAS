#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

float rsasum(const int N, const float* X, const int incX) {
  float_indexed *asumi = sialloc(SIDEFAULTFOLD);
  float asum;

  sisetzero(SIDEFAULTFOLD, asumi);

  sisasum(SIDEFAULTFOLD, N, X, incX, asumi);

  asum = ssiconv(SIDEFAULTFOLD, asumi);
  free(asumi);
  return asum;
}

