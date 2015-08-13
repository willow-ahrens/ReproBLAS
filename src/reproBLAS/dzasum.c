#include <reproBLAS.h>

#include "../../config.h"

double reproBLAS_dzasum(const int N, const void* X, const int incX) {
  return reproBLAS_rdzasum(DIDEFAULTFOLD, N, X, incX);
}
