#include <reproBLAS.h>

#include "../../config.h"

double reproBLAS_dasum(const int N, const double* X, const int incX) {
  return reproBLAS_rdasum(DIDEFAULTFOLD, N, X, incX);
}
