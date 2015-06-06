#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

void rcsum_sub(const int N, const void* X, const int incX, void *sum) {
  float_complex_indexed *sumi = cialloc(DEFAULT_FOLD);

  cisetzero(DEFAULT_FOLD, sumi);

  cicsum(DEFAULT_FOLD, N, X, incX, sumi);

  cciconv_sub(DEFAULT_FOLD, sumi, sum);
  free(sumi);
  return;
}

