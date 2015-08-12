#include <reproBLAS.h>
#include <idxdBLAS.h>

#include "../../config.h"

void rzdotu_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotu) {
  double_complex_indexed *dotui = idxd_zialloc(DIDEFAULTFOLD);

  idxd_zisetzero(DIDEFAULTFOLD, dotui);

  idxdBLAS_zizdotu(DIDEFAULTFOLD, N, X, incX, Y, incY, dotui);

  idxd_zziconv_sub(DIDEFAULTFOLD, dotui, dotu);
  free(dotui);
  return;
}

