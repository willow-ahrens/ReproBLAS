#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

float rssum(const int N, const float* X, const int incX) {
  float_indexed *sumi = sialloc(SIDEFAULTFOLD);
  float sum;

  sisetzero(SIDEFAULTFOLD, sumi);

  sissum(SIDEFAULTFOLD, N, X, incX, sumi);

  sum = ssiconv(SIDEFAULTFOLD, sumi);
  free(sumi);
  return sum;
}

