/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"

//TODO adjust this to reflect MAX_FOLD
#define BOUNDS_SIZE      64
#define BOUND_ZERO_INDEX 32
#define BIN_WIDTH        40
#define PREC             53

static double bounds[BOUNDS_SIZE];                   //initialized in bounds_initialize
static int    bounds_initialized = 0;                //initialized in bounds_initialize
static int    bound_min_index    = BOUND_ZERO_INDEX; //initialized in bounds_initialize
static int    bound_max_index    = BOUND_ZERO_INDEX; //initialized in bounds_initialize

int diwidth() {
  return BIN_WIDTH;
}

int dicapacity() {
  return 1 << (PREC - BIN_WIDTH - 2);
}

static void bounds_initialize() {
  int exp;
  int index;
  double step;

  if (bounds_initialized) {
    return;
  }

  bounds[BOUND_ZERO_INDEX] = 1.5;
  step = ldexp(1, BIN_WIDTH);

  exp = -1;
  index = BOUND_ZERO_INDEX + 1;
  while (exp * BIN_WIDTH  >= DBL_MIN_EXP) {
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
  while (exp * BIN_WIDTH <= DBL_MAX_EXP) {
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

int dmindex(double *repX){
  int index;

  bounds_initialize();

  if(isinf(repX[0])){
    index = bound_min_index;
  } else if(repX[0] == 0){
    index = bound_max_index;
  } else {
    frexp(repX[0], &index);
    index--;
    index /= BIN_WIDTH;
    index = BOUND_ZERO_INDEX - index;
  }
  return index;
}

int diindex(double_indexed *X){
  return dmindex(X);
}

int dindex(double X){
  int index;

  bounds_initialize();

  if(isinf(X)){
    index = bound_min_index;
  }else if(X == 0){
    index = bound_max_index;
  }else{
    frexp(X, &index);
    index += PREC - BIN_WIDTH;
    if(index < 0){
      index -= BIN_WIDTH - 1; //we want to round towards -infinity
    }
    index /= BIN_WIDTH;
    index = BOUND_ZERO_INDEX - 1 - index;
  }
  return index;
}

double dbound(int index){
  bounds_initialize();

  return bounds[index];
}

void dmbound(const int fold, int index, double *repY, int increpY) {
  int i;

  bounds_initialize();

  for (i = 0; i < fold; i++) {
    repY[i * increpY] = bounds[index + i];
  }
}
