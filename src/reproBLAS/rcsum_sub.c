#include <reproBLAS.h>
#include <idxdBLAS.h>

void reproBLAS_rcsum_sub(const int fold, const int N, const void* X, const int incX, void *sum) {
  float_complex_indexed *sumi = idxd_cialloc(fold);

  idxd_cisetzero(fold, sumi);

  idxdBLAS_cicsum(fold, N, X, incX, sumi);

  idxd_cciconv_sub(fold, sumi, sum);
  free(sumi);
  return;
}
