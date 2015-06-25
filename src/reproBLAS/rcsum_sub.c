#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

void rcsum_sub(const int N, const void* X, const int incX, void *sum) {
  float_complex_indexed *sumi = cialloc(SIDEFAULTFOLD);

  cisetzero(SIDEFAULTFOLD, sumi);

  cicsum(SIDEFAULTFOLD, N, X, incX, sumi);

  cciconv_sub(SIDEFAULTFOLD, sumi, sum);
  free(sumi);
  return;
}

