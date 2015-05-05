/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

double rdzasum(const int N, const void* X, const int incX) {
  double_indexed *asumi = dialloc(DEFAULT_FOLD);
  double asum;

  disetzero(DEFAULT_FOLD, asumi);

  dizasum(DEFAULT_FOLD, N, X, incX, asumi);

  asum = ddiconv(DEFAULT_FOLD, asumi);
  free(asumi);
  return asum;
}

