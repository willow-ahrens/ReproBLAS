#include <math.h>

#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

float rsnrm2(const int N, const float* X, const int incX) {
  float_indexed *ssq = idxd_sialloc(SIDEFAULTFOLD);
  float scl;
  float nrm2;

  idxd_sisetzero(SIDEFAULTFOLD, ssq);

  scl = sisssq(SIDEFAULTFOLD, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(idxd_ssiconv(SIDEFAULTFOLD, ssq));
  free(ssq);
  return nrm2;
}
