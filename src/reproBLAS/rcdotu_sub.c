#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

void rcdotu_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotu) {
  float_complex_indexed *dotui = cialloc(DEFAULT_FOLD);

  cisetzero(DEFAULT_FOLD, dotui);

  cicdotu(DEFAULT_FOLD, N, X, incX, Y, incY, dotui);

  cciconv_sub(DEFAULT_FOLD, dotui, dotu);
  free(dotui);
  return;
}

