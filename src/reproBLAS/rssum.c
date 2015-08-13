#include <reproBLAS.h>
#include <idxdBLAS.h>

float reproBLAS_rssum(const int fold, const int N, const float* X, const int incX) {
  float_indexed *sumi = idxd_sialloc(fold);
  float sum;

  idxd_sisetzero(fold, sumi);

  idxdBLAS_sissum(fold, N, X, incX, sumi);

  sum = idxd_ssiconv(fold, sumi);
  free(sumi);
  return sum;
}
