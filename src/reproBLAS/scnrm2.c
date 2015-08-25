#include <reproBLAS.h>

#include "../../config.h"

float reproBLAS_scnrm2(const int N, const void* X, const int incX) {
  return reproBLAS_rscnrm2(SIDEFAULTFOLD, N, X, incX);
}
