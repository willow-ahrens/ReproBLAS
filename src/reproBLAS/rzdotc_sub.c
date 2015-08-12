#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

void rzdotc_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  double_complex_indexed *dotci = idxd_zialloc(DIDEFAULTFOLD);

  idxd_zisetzero(DIDEFAULTFOLD, dotci);

  idxdBLAS_zizdotc(DIDEFAULTFOLD, N, X, incX, Y, incY, dotci);

  idxd_zziconv_sub(DIDEFAULTFOLD, dotci, dotc);
  free(dotci);
  return;
}

