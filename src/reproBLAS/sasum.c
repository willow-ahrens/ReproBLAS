#include <reproBLAS.h>

#include "../../config.h"

float reproBLAS_sasum(const int N, const float* X, const int incX) {
  return reproBLAS_rsasum(SIDEFAULTFOLD, N, X, incX);
}
