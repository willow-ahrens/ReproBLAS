#include <reproBLAS.h>
#include <idxdBLAS.h>

void reproBLAS_rzdotu_sub(const int fold, const int N, const void* X, const int incX, const void *Y, const int incY, void *dotu) {
  double_complex_indexed *dotui = idxd_zialloc(fold);

  idxd_zisetzero(fold, dotui);

  idxdBLAS_zizdotu(fold, N, X, incX, Y, incY, dotui);

  idxd_zziconv_sub(fold, dotui, dotu);
  free(dotui);
  return;
}
