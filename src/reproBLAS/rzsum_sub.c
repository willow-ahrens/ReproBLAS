#include <reproBLAS.h>
#include <idxdBLAS.h>

#include "../../config.h"

void reproBLAS_rzsum_sub(const int N, const void* X, const int incX, void *sum) {
  double_complex_indexed *sumi = idxd_zialloc(DIDEFAULTFOLD);

  idxd_zisetzero(DIDEFAULTFOLD, sumi);

  idxdBLAS_zizsum(DIDEFAULTFOLD, N, X, incX, sumi);

  idxd_zziconv_sub(DIDEFAULTFOLD, sumi, sum);
  free(sumi);
  return;
}

