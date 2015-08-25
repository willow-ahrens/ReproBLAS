#include <reproBLAS.h>
#include <idxdBLAS.h>

void reproBLAS_rcdotu_sub(const int fold, const int N, const void* X, const int incX, const void *Y, const int incY, void *dotu) {
  float_complex_indexed *dotui = idxd_cialloc(fold);

  idxd_cisetzero(fold, dotui);

  idxdBLAS_cicdotu(fold, N, X, incX, Y, incY, dotui);

  idxd_cciconv_sub(fold, dotui, dotu);
  free(dotui);
  return;
}
