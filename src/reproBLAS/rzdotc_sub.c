#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

void rzdotc_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  double_complex_indexed *dotci = zialloc(DEFAULT_FOLD);

  zisetzero(DEFAULT_FOLD, dotci);

  zizdotc(DEFAULT_FOLD, N, X, incX, Y, incY, dotci);

  zziconv_sub(DEFAULT_FOLD, dotci, dotc);
  free(dotci);
  return;
}

