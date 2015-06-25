#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

void rzsum_sub(const int N, const void* X, const int incX, void *sum) {
  double_complex_indexed *sumi = zialloc(DIDEFAULTFOLD);

  zisetzero(DIDEFAULTFOLD, sumi);

  zizsum(DIDEFAULTFOLD, N, X, incX, sumi);

  zziconv_sub(DIDEFAULTFOLD, sumi, sum);
  free(sumi);
  return;
}

