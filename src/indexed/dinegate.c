#include "indexed.h"

void dmnegate(const int fold, double* repX, const int increpX, double* carX, const int inccarX) {
  int i;
  for (i = 0; i < fold; i++) {
    repX[i * increpX] = (3 * ufp(repX[i * increpX])) - repX[i * increpX];
    carX[i * inccarX] = -carX[i * inccarX];
  }
}

void dinegate(const int fold, double_indexed* X){
  dmnegate(fold, X, 1, X + fold, 1);
}

void zmnegate(const int fold, double* repX, const int increpX, double* carX, const int inccarX) {
  dmnegate(fold, repX, 2 * increpX, carX, 2 * inccarX);
  dmnegate(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX);
}

void zinegate(const int fold, double_complex_indexed* X){
  zmnegate(fold, X, 1, X + 2 * fold, 1);
}
