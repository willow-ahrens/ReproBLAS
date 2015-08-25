#include <reproBLAS.h>
#include <idxdBLAS.h>

void reproBLAS_rzsum_sub(const int fold, const int N, const void* X, const int incX, void *sum) {
  double_complex_indexed *sumi = idxd_zialloc(fold);

  idxd_zisetzero(fold, sumi);

  idxdBLAS_zizsum(fold, N, X, incX, sumi);

  idxd_zziconv_sub(fold, sumi, sum);
  free(sumi);
  return;
}
