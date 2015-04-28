#include <stdio.h>

#include "indexed.h"

void zmprint(const int fold, double *repX, int increpX, double *carX, int inccarX){
  int i;
  double M;
  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    M = ufp(repX[0]);
    printf("(2^%d: %g #%g =%g", (int)log2(M), repX[0] - 1.5*M, carX[0], (carX[0] - 6) * 0.25 * M + repX[0]);
    M = ufp(repX[1]);
    printf("| 2^%d: %g #%g =%g)\n", (int)log2(M), repX[1] - 1.5*M, carX[1], (carX[1] - 6) * 0.25 * M + repX[1]);
  }
}

void ziprint(const int fold, double_complex_indexed *X){
  zmprint(fold, X, 2, X + 2 * fold, 2);
}
