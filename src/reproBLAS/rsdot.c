#include <reproBLAS.h>
#include <idxdBLAS.h>

#include "../../config.h"

float reproBLAS_rsdot(const int N, const float* X, const int incX, const float *Y, const int incY) {
  float_indexed *doti = idxd_sialloc(SIDEFAULTFOLD);
  float dot;

  idxd_sisetzero(SIDEFAULTFOLD, doti);

  idxdBLAS_sisdot(SIDEFAULTFOLD, N, X, incX, Y, incY, doti);

  dot = idxd_ssiconv(SIDEFAULTFOLD, doti);
  free(doti);
  return dot;
}

