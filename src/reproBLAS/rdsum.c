#include <reproBLAS.h>
#include <idxdBLAS.h>

double reproBLAS_rdsum(const int fold, const int N, const double* X, const int incX) {
  double_indexed *sumi = idxd_dialloc(fold);
  double sum;

  idxd_disetzero(fold, sumi);

  idxdBLAS_didsum(fold, N, X, incX, sumi);

  sum = idxd_ddiconv(fold, sumi);
  free(sumi);
  return sum;
}
