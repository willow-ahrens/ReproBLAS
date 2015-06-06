#include <math.h>

#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

float rscnrm2(const int N, const void* X, const int incX) {
  float_indexed *ssq = sialloc(DEFAULT_FOLD);
  float scl;
  float nrm2;

  sisetzero(DEFAULT_FOLD, ssq);

  scl = sicssq(DEFAULT_FOLD, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(ssiconv(DEFAULT_FOLD, ssq));
  free(ssq);
  return nrm2;
}
