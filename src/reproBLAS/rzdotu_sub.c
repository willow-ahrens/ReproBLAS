#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

void rzdotu_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotu) {
  double_complex_indexed *dotui = zialloc(DIDEFAULTFOLD);

  zisetzero(DIDEFAULTFOLD, dotui);

  zizdotu(DIDEFAULTFOLD, N, X, incX, Y, incY, dotui);

  zziconv_sub(DIDEFAULTFOLD, dotui, dotu);
  free(dotui);
  return;
}

