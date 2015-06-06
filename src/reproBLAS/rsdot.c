#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

float rsdot(const int N, const float* X, const int incX, const float *Y, const int incY) {
  float_indexed *doti = sialloc(DEFAULT_FOLD);
  float dot;

  sisetzero(DEFAULT_FOLD, doti);

  sisdot(DEFAULT_FOLD, N, X, incX, Y, incY, doti);

  dot = ssiconv(DEFAULT_FOLD, doti);
  free(doti);
  return dot;
}

