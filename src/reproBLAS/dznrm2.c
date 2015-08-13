#include <reproBLAS.h>

#include "../../config.h"

double reproBLAS_dznrm2(const int N, const void* X, const int incX) {
  return reproBLAS_rdznrm2(DIDEFAULTFOLD, N, X, incX);
}
