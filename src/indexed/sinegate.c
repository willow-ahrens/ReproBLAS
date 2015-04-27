#include "indexed.h"

void smnegate(const int fold, float* repX, const int increpX, float* carX, const int inccarX) {
  int i;
  for (i = 0; i < fold; i++) {
    repX[i * increpX] = (3 * ufp(repX[i * increpX])) - repX[i * increpX];
    carX[i * inccarX] = -carX[i * inccarX];
  }
}

void sinegate(const int fold, float_indexed* X){
  smnegate(fold, X, 1, X + fold, 1);
}

void cmnegate(const int fold, float* repX, const int increpX, float* carX, const int inccarX) {
  smnegate(fold, repX, 2 * increpX, carX, 2 * inccarX);
  smnegate(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX);
}

void cinegate(const int fold, float_complex_indexed* X){
  cmnegate(fold, X, 1, X + 2 * fold, 1);
}
