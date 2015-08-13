#include <math.h>

#include <reproBLAS.h>
#include <idxdBLAS.h>

float reproBLAS_rscnrm2(const int fold, const int N, const void* X, const int incX) {
  float_indexed *ssq = idxd_sialloc(fold);
  float scl;
  float nrm2;

  idxd_sisetzero(fold, ssq);

  scl = idxdBLAS_sicssq(fold, N, X, incX, 0.0, ssq);

  nrm2 = scl * sqrt(idxd_ssiconv(fold, ssq));
  free(ssq);
  return nrm2;
}
