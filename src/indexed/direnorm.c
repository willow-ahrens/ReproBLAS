/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

void dmrenorm(double* repX, int increpX, double* carX, int inccarX, int fold) {
  int i;
  double M;
  double repX0;
  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    repX0 = repX[0];
    if (repX0 == 0.0)
      continue;

    M = ufp(repX0);
    if (repX0 >= (M * 1.75)) {
      repX[0] -= M * 0.25;
      carX[0] += 1;
    }
    else if (repX0 < (M * 1.25)) {
      repX[0] += M * 0.5;
      carX[0] -= 2;
    }
    else if (repX0 < (M * 1.5)) {
      repX[0] += M * 0.25;
      carX[0] -= 1;
    }
  }
}

void direnorm(double_indexed *X, int fold) {
  dmrenorm(X, 1, X + fold, 1, fold);
}

void zmrenorm(double* repX, int increpX, double* carX, int inccarX, int fold) {
  dmrenorm(repX, 2 * increpX, carX, 2 * inccarX, fold);
  dmrenorm(repX + 1, 2 * increpX, carX + 1, 2 * inccarX, fold);
}

void zirenorm(double_complex_indexed *X, int fold) {
  zmrenorm(X, 1, X + 2 * fold, 1, fold);
}
