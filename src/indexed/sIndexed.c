/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"

#define BOUNDS_SIZE      20
#define BOUND_ZERO_INDEX 10
#define BIN_WIDTH        15
#define PREC             23

static float bounds[BOUNDS_SIZE];     //initialized in bounds_initialize
static int   bounds_initialized = 0;  //initialized in bounds_initialize
static int   bound_min_index    = 32; //initialized in bounds_initialize
static int   bound_max_index    = 32; //initialized in bounds_initialize

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

int siindex(float_indexed *X){
  int index;

  bounds_initialize();

  if(isinf(X[0])){
    index = bound_min_index;
  } else if(X[0] == 0){
    index = bound_max_index;
  } else {
    frexpf(X[0], &index);
    index--;
    index /= BIN_WIDTH;
    index = BOUND_ZERO_INDEX - index;
  }
  return index;
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
    index += PREC - BIN_WIDTH - 1;
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

void sIprint1(int n, float* x, float* carry, int inc) {
	int i;
	float M;
	for (i = 0; i < n; i++) {
		M = ufpf(x[i*inc]);
		printf("{M:2^%g # %g :: %g (%g)}", log2f(M), carry[i*inc], x[i*inc] - 1.5*M, x[i*inc] + (carry[i*inc]-6)*0.25*M);
	}
}

void cIprint1(int n, float complex* x, float* carry, int inc) {
	/*
	printf("\nReal: ");
	sIprint_(n, (float*) x, carry, inc*2);
	printf("\nImaginary: ");
	sIprint_(n, ((float*) x) + 1, carry+1, inc*2);
	*/
	int i;
	float M;
	float* ptr = (float*) x;
	for (i = 0; i < n; i++, ptr += 2 * inc, carry += 2 * inc) {
		M = ufpf(ptr[0]);
		printf("M:2^%2g", log2f(M));
		printf("# %4g", carry[0]);
		printf(" %8.3g (%8.3g) ", ptr[0], ptr[0] - 1.5 * M);

		M = ufpf(ptr[1]);
		printf(" || M:2^%2g", log2f(M));
		printf("# %6g", carry[1]);
		printf(" %8.3g (%8.3g) ", ptr[1], ptr[1] - 1.5 * M);
		printf("\n");
	}
}
