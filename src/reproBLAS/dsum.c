#include <reproBLAS.h>

#include "../../config.h"

double reproBLAS_dsum(const int N, const double* X, const int incX) {
  return reproBLAS_rdsum(DIDEFAULTFOLD, N, X, incX);
}
