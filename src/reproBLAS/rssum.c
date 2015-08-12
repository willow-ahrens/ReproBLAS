#include <reproBLAS.h>
#include <idxdBLAS.h>

#include "../../config.h"

float rssum(const int N, const float* X, const int incX) {
  float_indexed *sumi = idxd_sialloc(SIDEFAULTFOLD);
  float sum;

  idxd_sisetzero(SIDEFAULTFOLD, sumi);

  idxdBLAS_sissum(SIDEFAULTFOLD, N, X, incX, sumi);

  sum = idxd_ssiconv(SIDEFAULTFOLD, sumi);
  free(sumi);
  return sum;
}

