#include <stdio.h>

#include "indexed.h"

void dmprint(double *repX, int increpX, double *carX, int inccarX, int fold) {
  int i;
  double M;
  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    M = ufp(repX[0]);
    printf("(2^%d: %g #%g =%g)", (int)log2(M), repX[0] - 1.5*M, carX[0], (carX[0] - 6) * 0.25 * M + repX[0]);
  }
}

void diprint(double_indexed *X, int fold) {
  dmprint(X, 1, X + fold, 1, fold);
}
