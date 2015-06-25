#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

float rsdot(const int N, const float* X, const int incX, const float *Y, const int incY) {
  float_indexed *doti = sialloc(SIDEFAULTFOLD);
  float dot;

  sisetzero(SIDEFAULTFOLD, doti);

  sisdot(SIDEFAULTFOLD, N, X, incX, Y, incY, doti);

  dot = ssiconv(SIDEFAULTFOLD, doti);
  free(doti);
  return dot;
}

