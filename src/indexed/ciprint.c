#include <stdio.h>

#include "indexed.h"

void cmprint(float* repX, int increpX, float* carX, int inccarX, int fold) {
  int i;
  float M;
  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    M = ufpf(repX[0]);
    printf("(%g)", repX[0]);
    printf("(2^%d: %g #%g =%g", (int)log2f(M), repX[0] - 1.5*M, carX[0], (carX[0] - 6) * 0.25 * M + repX[0]);
    M = ufpf(repX[1]);
    printf("(%g)", repX[1]);
    printf("| 2^%d: %g #%g =%g)\n", (int)log2f(M), repX[1] - 1.5*M, carX[1], (carX[1] - 6) * 0.25 * M + repX[1]);
  }
}

void ciprint(float_complex_indexed *X, int fold) {
  cmprint(X, 2, X + fold * 2, 2, fold);
}
