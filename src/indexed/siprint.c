#include <stdio.h>

#include "indexed.h"

void smprint(float *repX, int increpX, float *carX, int inccarX, int fold) {
  int i;
  float M;
  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    M = ufpf(repX[0]);
    printf("(2^%d: %g #%g =%g)\n", (int)log2f(M), repX[0] - 1.5*M, carX[0], (carX[0] - 6) * 0.25 * M + repX[0]);
  }
}

void siprint(float_indexed *X, int fold) {
  smprint(X, 1, X + fold, 1, fold);
}
