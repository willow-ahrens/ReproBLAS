#include <reproBLAS.h>
#include <idxdBLAS.h>

double reproBLAS_rddot(const int fold, const int N, const double* X, const int incX, const double *Y, const int incY) {
  double_indexed *doti = idxd_dialloc(fold);
  double dot;

  idxd_disetzero(fold, doti);

  idxdBLAS_diddot(fold, N, X, incX, Y, incY, doti);

  dot = idxd_ddiconv(fold, doti);
  free(doti);
  return dot;
}
