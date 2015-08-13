#include <reproBLAS.h>

#include "../../config.h"

float reproBLAS_sdot(const int N, const float* X, const int incX, const float *Y, const int incY) {
  return reproBLAS_rsdot(SIDEFAULTFOLD, N, X, incX, Y, incY);
}
