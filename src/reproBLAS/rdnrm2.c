#include <math.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

double reproBLAS_rdnrm2(const int fold, const int N, const double* X, const int incX) {
  double_indexed *ssq = idxd_dialloc(fold);
  double scl;
  double nrm2;

  idxd_disetzero(fold, ssq);

  scl = idxdBLAS_didssq(fold, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(idxd_ddiconv(fold, ssq));
  free(ssq);
  return nrm2;
}
