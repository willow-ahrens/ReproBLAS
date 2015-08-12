#include <reproBLAS.h>
#include <idxdBLAS.h>

#include "../../config.h"

void rcdotc_sub(const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  float_complex_indexed *dotci = idxd_cialloc(SIDEFAULTFOLD);

  idxd_cisetzero(SIDEFAULTFOLD, dotci);

  idxdBLAS_cicdotc(SIDEFAULTFOLD, N, X, incX, Y, incY, dotci);

  idxd_cciconv_sub(SIDEFAULTFOLD, dotci, dotc);
  free(dotci);
  return;
}

