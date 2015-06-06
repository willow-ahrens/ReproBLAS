#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

float rscasum(const int N, const void* X, const int incX) {
  float_indexed *asumi = sialloc(DEFAULT_FOLD);
  float asum;

  sisetzero(DEFAULT_FOLD, asumi);

  sicasum(DEFAULT_FOLD, N, X, incX, asumi);

  asum = ssiconv(DEFAULT_FOLD, asumi);
  free(asumi);
  return asum;
}

