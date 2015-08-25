#include <reproBLAS.h>

#include "../../config.h"

float reproBLAS_snrm2(const int N, const float* X, const int incX) {
  return reproBLAS_rsnrm2(SIDEFAULTFOLD, N, X, incX);
}
