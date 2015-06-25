#include <math.h>

#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

float rscnrm2(const int N, const void* X, const int incX) {
  float_indexed *ssq = sialloc(SIDEFAULTFOLD);
  float scl;
  float nrm2;

  sisetzero(SIDEFAULTFOLD, ssq);

  scl = sicssq(SIDEFAULTFOLD, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(ssiconv(SIDEFAULTFOLD, ssq));
  free(ssq);
  return nrm2;
}
