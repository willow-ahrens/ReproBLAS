#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

float rsasum(const int N, const float* X, const int incX) {
  float_indexed *asumi = sialloc(DEFAULT_FOLD);
  float asum;

  sisetzero(DEFAULT_FOLD, asumi);

  sisasum(DEFAULT_FOLD, N, X, incX, asumi);

  asum = ssiconv(DEFAULT_FOLD, asumi);
  free(asumi);
  return asum;
}

