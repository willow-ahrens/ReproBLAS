#include <reproBLAS.h>
#include <idxdBLAS.h>

#include "../../config.h"

void reproBLAS_rcsum_sub(const int N, const void* X, const int incX, void *sum) {
  float_complex_indexed *sumi = idxd_cialloc(SIDEFAULTFOLD);

  idxd_cisetzero(SIDEFAULTFOLD, sumi);

  idxdBLAS_cicsum(SIDEFAULTFOLD, N, X, incX, sumi);

  idxd_cciconv_sub(SIDEFAULTFOLD, sumi, sum);
  free(sumi);
  return;
}

