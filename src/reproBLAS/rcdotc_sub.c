#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

void rcdotc_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  float_complex_indexed *dotci = cialloc(SIDEFAULTFOLD);

  cisetzero(SIDEFAULTFOLD, dotci);

  cicdotc(SIDEFAULTFOLD, N, X, incX, Y, incY, dotci);

  cciconv_sub(SIDEFAULTFOLD, dotci, dotc);
  free(dotci);
  return;
}

