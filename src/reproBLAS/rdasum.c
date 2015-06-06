#include <reproBLAS.h>
#include <indexedBLAS.h>

#include "../../config.h"

double rdasum(const int N, const double* X, const int incX) {
  double_indexed *asumi = dialloc(DEFAULT_FOLD);
  double asum;

  disetzero(DEFAULT_FOLD, asumi);

  didasum(DEFAULT_FOLD, N, X, incX, asumi);

  asum = ddiconv(DEFAULT_FOLD, asumi);
  free(asumi);
  return asum;
}

