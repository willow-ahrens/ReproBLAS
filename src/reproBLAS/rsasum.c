#include <reproBLAS.h>
#include <idxdBLAS.h>

float reproBLAS_rsasum(const int fold, const int N, const float* X, const int incX) {
  float_indexed *asumi = idxd_sialloc(fold);
  float asum;

  idxd_sisetzero(fold, asumi);

  idxdBLAS_sisasum(fold, N, X, incX, asumi);

  asum = idxd_ssiconv(fold, asumi);
  free(asumi);
  return asum;
}
