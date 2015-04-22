/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"

#define BOUNDS_SIZE      64
#define BOUND_ZERO_INDEX 32
#define BIN_WIDTH        40
#define PREC             53

static double bounds[BOUNDS_SIZE];     //initialized in bounds_initialize
static int    bounds_initialized = 0;  //initialized in bounds_initialize
static int    bound_min_index    = 32; //initialized in bounds_initialize
static int    bound_max_index    = 32; //initialized in bounds_initialize

int dIWidth() {
  return BIN_WIDTH;
}

int dICapacity() {
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

double dbound(int index){
  bounds_initialize();

  return bounds[index];
}

// COMPUTE THE BOUNDARIES BASED ON MAXIMUM ABSOLUTE VALUE
int dIBoundary(int fold, double max, double* M, int inc) {
  int i;
  int index;

  bounds_initialize();

  index = dindex(max);

  for (i = 0; i < fold; i++) {
    M[i * inc] = bounds[index + i];
  }
  return index;
}

int diindex(double_indexed *x){
  int index;

  bounds_initialize();

  if(isinf(x[0])){
    index = bound_min_index;
  } else if(x[0] == 0){
    index = bound_max_index;
  } else {
    frexp(x[0], &index);
    index--;
    index /= BIN_WIDTH;
    index = BOUND_ZERO_INDEX - index;
  }
  return index;
}

int dindex(double x){
  int index;

  bounds_initialize();

  if(isinf(x)){
    index = bound_min_index;
  } else if(x == 0){
    index = bound_max_index;
  } else {
    frexp(x, &index);
    index += PREC - BIN_WIDTH - 1;
    if(index < 0){
      index -= BIN_WIDTH - 1; //we want to round towards -infinity
    }
    index /= BIN_WIDTH;
    index = BOUND_ZERO_INDEX - 1 - index;
  }
  return index;
}



void dIprint1(int n, double *x, double *c, int inc) {
int i;
double M;
for (i = 0; i < n; i++, x += inc, c += inc) {
M = ufp(x[0]);
printf("{M:2^%g # %g :: %g (%.16g)}", log2(M), c[0], x[0] - 1.5*M,
(c[0] - 6) * 0.25 * M + x[0]);
}
}

void zIprint1(int n, double complex* x, double complex* carry, int inc) {
int i;
double M;
double* ptr = (double*) x;
double* cptr = (double*) carry;
inc *= 2;
for (i = 0; i < n; i++, ptr += inc, cptr += inc) {
M = ufp(ptr[0]);
printf("M:2^%2g", log2(M));
printf("# %4g", cptr[0]);
printf(" %.8g (%8.3g) ", ptr[0], ptr[0] - 1.5 * M);

M = ufp(ptr[1]);
printf(" || M:2^%2g", log2(M));
printf("# %6g", cptr[1]);
printf(" %.8g (%8.3g) ", ptr[1], ptr[1] - 1.5 * M);
printf("\n");
}
}
