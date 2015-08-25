#include <reproBLAS.h>

#include "../../config.h"

double reproBLAS_ddot(const int N, const double* X, const int incX, const double *Y, const int incY) {
  return reproBLAS_rddot(DIDEFAULTFOLD, N, X, incX, Y, incY);
}
