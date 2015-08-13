#include <reproBLAS.h>
#include <idxdBLAS.h>

double reproBLAS_rdasum(const int fold, const int N, const double* X, const int incX) {
  double_indexed *asumi = idxd_dialloc(fold);
  double asum;

  idxd_disetzero(fold, asumi);

  idxdBLAS_didasum(fold, N, X, incX, asumi);

  asum = idxd_ddiconv(fold, asumi);
  free(asumi);
  return asum;
}
