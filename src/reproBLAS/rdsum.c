#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

double rdsum(const int N, const double* X, const int incX) {
  double_indexed *sumi = dialloc(DEFAULT_FOLD);
  double sum;

  disetzero(DEFAULT_FOLD, sumi);

  didsum(DEFAULT_FOLD, N, X, incX, sumi);

  sum = ddiconv(DEFAULT_FOLD, sumi);
  free(sumi);
  return sum;
}

