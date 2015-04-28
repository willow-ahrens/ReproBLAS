#include <stdlib.h>
#include <string.h>

#include "indexed.h"

size_t sisize(const int fold){
  return 2*fold*sizeof(float);
}

size_t cisize(const int fold){
  return 4*fold*sizeof(float);
}

int sinum(const int fold){
  return 2*fold;
}

int cinum(const int fold){
  return 4*fold;
}

float_indexed *sialloc(const int fold){
  return (float_indexed*)malloc(sisize(fold));
}

float_complex_indexed *cialloc(const int fold){
  return (float_complex_indexed*)malloc(cisize(fold));
}

void smsmset(const int fold, const float *repX, const int increpX, const float *carX, const int inccarX, float *repY, const int increpY, float *carY, const int inccarY){
  int i;
  for(i = 0; i < fold; i++){
    repY[i * increpY] = repX[i * increpX];
    carY[i * inccarY] = carX[i * inccarX];
  }
}

void sisiset(const int fold, const float_indexed *X, float_indexed *Y){
  memcpy(Y, X, sisize(fold));
}

void cmcmset(const int fold, const float *repX, const int increpX, const float *carX, const int inccarX, float *repY, const int increpY, float *carY, const int inccarY){
  smsmset(fold, repX, 2 * increpX, carX, 2 * inccarX, repY, 2 * increpY, carY, 2 * inccarY);
  smsmset(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX, repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}

void ciciset(const int fold, const float_complex_indexed *X, float_complex_indexed *Y){
  memcpy(Y, X, cisize(fold));
}

void cmsmset(const int fold, const float *repX, const int increpX, const float *carX, const int inccarX, float *repY, const int increpY, float *carY, const int inccarY){
  smsmset(fold, repX, increpX, carX, inccarX, repY, 2 * increpY, carY, 2 * inccarY);
  smsetzero(fold, repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}

void cisiset(const int fold, const float_indexed *X, float_complex_indexed *Y){
  cmsmset(fold, X, 1, X + fold, 1, Y, 1, Y + 2 * fold, 1);
}

void smsetzero(const int fold, float *repX, const int increpX, float *carX, const int inccarX){
  int i;
  for(i = 0; i < fold; i++){
    repX[i * increpX] = 0.0;
    carX[i * inccarX] = 0.0;
  }
}

void sisetzero(const int fold, float_indexed *X){
  memset(X, 0, sisize(fold));
}

void cmsetzero(const int fold, float *repX, const int increpX, float *carX, const int inccarX){
  smsetzero(fold, repX, 2 * increpX, carX, 2 * inccarX);
  smsetzero(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX);
}

void cisetzero(const int fold, float_complex_indexed *X){
  memset(X, 0, cisize(fold));
}
