#include <reproBLAS.h>
#include <idxdBLAS.h>

void reproBLAS_rcdotc_sub(const int fold, const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  float_complex_indexed *dotci = idxd_cialloc(fold);

  idxd_cisetzero(fold, dotci);

  idxdBLAS_cicdotc(fold, N, X, incX, Y, incY, dotci);

  idxd_cciconv_sub(fold, dotci, dotc);
  free(dotci);
  return;
}
