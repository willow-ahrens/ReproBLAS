/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"

#define BOUNDS_SIZE      32
#define BOUND_ZERO_INDEX 16
#define BIN_WIDTH        13
#define PREC             24

static float bounds[BOUNDS_SIZE];                   //initialized in bounds_initialize
static int   bounds_initialized = 0;                //initialized in bounds_initialize
static int   bound_min_index    = BOUND_ZERO_INDEX; //initialized in bounds_initialize
static int   bound_max_index    = BOUND_ZERO_INDEX; //initialized in bounds_initialize

int sIWidth() {
  return BIN_WIDTH;
}

int sICapacity() {
  return 1 << (PREC - BIN_WIDTH - 2);
}

static void bounds_initialize() {
  int exp;
  int index;
  float step;

  if (bounds_initialized) {
    return;
  }

  bounds[BOUND_ZERO_INDEX] = 1.5;
  step = ldexpf(1, BIN_WIDTH);

  exp = -1;
  index = BOUND_ZERO_INDEX + 1;
  while (exp * BIN_WIDTH  >= FLT_MIN_EXP) {
    bounds[index] = bounds[index - 1] / step;
    index++;
    exp--;
  }
  bound_max_index = index;
  while (index < BOUNDS_SIZE) {
    bounds[index] = 0.0;
    index++;
  }

  exp = 1;
  index = BOUND_ZERO_INDEX - 1;
  while (exp * BIN_WIDTH <= FLT_MAX_EXP) {
    bounds[index] = bounds[index + 1] * step;
    index--;
    exp++;
  }
  bound_min_index = index;
  while (index >= 0) {
    bounds[index] = bounds[bound_min_index + 1] * step;
    index--;
  }

  bounds_initialized = 1;
}

int smindex(float *repX){
  int index;

  bounds_initialize();

  if(isinf(repX[0])){
    index = bound_min_index;
  } else if(repX[0] == 0){
    index = bound_max_index;
  } else {
    frexpf(repX[0], &index);
    index--;
    index /= BIN_WIDTH;
    index = BOUND_ZERO_INDEX - index;
  }
  return index;
}

int siindex(float_indexed *X){
  return smindex(X);
}

int sindex(float X){
  int index;

  bounds_initialize();

  if(isinf(X)){
    index = bound_min_index;
  }else if(X == 0){
    index = bound_max_index;
  }else{
    frexpf(X, &index);
    index += PREC - BIN_WIDTH;
    if(index < 0){
      index -= BIN_WIDTH - 1; //we want to round towards -infinity
    }
    index /= BIN_WIDTH;
    index = BOUND_ZERO_INDEX - 1 - index;
  }
  return index;
}

double sbound(int index){
  bounds_initialize();

  return bounds[index];
}

void smbound(int index, float *repY, int increpY, int fold) {
  int i;

  bounds_initialize();

  for (i = 0; i < fold; i++) {
    repY[i * increpY] = bounds[index + i];
  }
}
