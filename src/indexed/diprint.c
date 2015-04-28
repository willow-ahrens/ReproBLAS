#include <stdio.h>

#include "indexed.h"

void dmprint(const int fold, double *repX, int increpX, double *carX, int inccarX) {
  int i;
  double M;
  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    M = ufp(repX[0]);
    printf("(2^%d: %g #%g =%g)", (int)log2(M), repX[0] - 1.5*M, carX[0], (carX[0] - 6) * 0.25 * M + repX[0]);
  }
}

void diprint(const int fold, double_indexed *X) {
  dmprint(fold, X, 1, X + fold, 1);
}
