#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

void rcdotc_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  float_complex_indexed *dotci = cialloc(DEFAULT_FOLD);

  cisetzero(DEFAULT_FOLD, dotci);

  cicdotc(DEFAULT_FOLD, N, X, incX, Y, incY, dotci);

  cciconv_sub(DEFAULT_FOLD, dotci, dotc);
  free(dotci);
  return;
}

