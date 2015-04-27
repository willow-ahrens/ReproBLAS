#include <stdlib.h>
#include <string.h>

#include "indexed.h"

size_t disize(const int fold){
  return 2*fold*sizeof(double);
}

size_t zisize(const int fold){
  return 4*fold*sizeof(double);
}

size_t dinum(const int fold){
  return 2*fold;
}

size_t zinum(const int fold){
  return 4*fold;
}

double_indexed *dialloc(const int fold){
  return (double_indexed*)malloc(disize(fold));
}

double_complex_indexed *zialloc(const int fold){
  return (double_complex_indexed*)malloc(zisize(fold));
}

void dmset(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double *repY, const int increpY, double *carY, const int inccarY){
  int i;
  for(i = 0; i < fold; i++){
    repY[i * increpY] = repX[i * increpX];
    carY[i * inccarY] = carX[i * inccarX];
  }
}

void diset(const int fold, const double_indexed *X, double_indexed *Y){
  memcpy(Y, X, disize(fold));
}

void zmset(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double *repY, const int increpY, double *carY, const int inccarY){
  dmset(fold, repX, 2 * increpX, carX, 2 * inccarX, repY, 2 * increpY, carY, 2 * inccarY);
  dmset(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX, repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}

void ziset(const int fold, const double_complex_indexed *X, double_complex_indexed *Y){
  memcpy(Y, X, zisize(fold));
}

void dmsetzero(const int fold, double *repX, const int increpX, double *carX, const int inccarX){
  int i;
  for(i = 0; i < fold; i++){
    repX[i * increpX] = 0.0;
    carX[i * inccarX] = 0.0;
  }
}

void disetzero(const int fold, double_indexed *X){
  memset(X, 0, disize(fold));
}

void zmsetzero(const int fold, double *repX, const int increpX, double *carX, const int inccarX){
  dmsetzero(fold, repX, 2 * increpX, carX, 2 * inccarX);
  dmsetzero(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX);
}

void zisetzero(const int fold, double_complex_indexed *X){
  memset(X, 0, zisize(fold));
}
