#include <reproBLAS.h>
#include <idxdBLAS.h>

void reproBLAS_rzdotc_sub(const int fold, const int N, const void* X, const int incX, const void *Y, const int incY, void *dotc) {
  double_complex_indexed *dotci = idxd_zialloc(fold);

  idxd_zisetzero(fold, dotci);

  idxdBLAS_zizdotc(fold, N, X, incX, Y, incY, dotci);

  idxd_zziconv_sub(fold, dotci, dotc);
  free(dotci);
  return;
}
