#include <reproBLAS.h>

#include "../../config.h"

double reproBLAS_dnrm2(const int N, const double* X, const int incX) {
  return reproBLAS_rdnrm2(DIDEFAULTFOLD, N, X, incX);
}
