#include <reproBLAS.h>

#include "../../config.h"

float reproBLAS_ssum(const int N, const float* X, const int incX) {
  return reproBLAS_rssum(SIDEFAULTFOLD, N, X, incX);
}
