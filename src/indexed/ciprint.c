#include <stdio.h>

#include "indexed.h"

void cmprint(const int fold, float* repX, int increpX, float* carX, int inccarX) {
  int i;
  float M;
  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    M = ufpf(repX[0]);
    printf("(2^%d: %g #%g =%g", (int)log2f(M), repX[0] - 1.5*M, carX[0], (carX[0] - 6) * 0.25 * M + repX[0]);
    M = ufpf(repX[1]);
    printf("| 2^%d: %g #%g =%g)\n", (int)log2f(M), repX[1] - 1.5*M, carX[1], (carX[1] - 6) * 0.25 * M + repX[1]);
  }
}

void ciprint(const int fold, float_complex_indexed *X) {
  cmprint(fold, X, 2, X + fold * 2, 2);
}
