#include <reproBLAS.h>

#include "../../config.h"

float reproBLAS_scasum(const int N, const void* X, const int incX) {
  return reproBLAS_rscasum(SIDEFAULTFOLD, N, X, incX);
}
