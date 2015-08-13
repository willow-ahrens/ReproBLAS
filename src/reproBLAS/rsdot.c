#include <reproBLAS.h>
#include <idxdBLAS.h>

float reproBLAS_rsdot(const int fold, const int N, const float* X, const int incX, const float *Y, const int incY) {
  float_indexed *doti = idxd_sialloc(fold);
  float dot;

  idxd_sisetzero(fold, doti);

  idxdBLAS_sisdot(fold, N, X, incX, Y, incY, doti);

  dot = idxd_ssiconv(fold, doti);
  free(doti);
  return dot;
}
