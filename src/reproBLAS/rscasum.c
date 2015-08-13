#include <reproBLAS.h>
#include <idxdBLAS.h>

float reproBLAS_rscasum(const int fold, const int N, const void* X, const int incX) {
  float_indexed *asumi = idxd_sialloc(fold);
  float asum;

  idxd_sisetzero(fold, asumi);

  idxdBLAS_sicasum(fold, N, X, incX, asumi);

  asum = idxd_ssiconv(fold, asumi);
  free(asumi);
  return asum;
}
