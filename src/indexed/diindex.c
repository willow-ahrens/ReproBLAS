/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../../config.h"

#define PREC             53
#define BIN_WIDTH        40
#define BOUNDS_SIZE      ((DBL_MAX_EXP - DBL_MIN_EXP)/BIN_WIDTH + MAX_FOLD + 2)
#define BOUNDS_ZERO_INDEX (DBL_MAX_EXP/BIN_WIDTH + 1)

static double bounds[BOUNDS_SIZE];                     //initialized in bounds_initialize
static int    bounds_initialized = 0;                  //initialized in bounds_initialize
static int    bounds_min_index    = BOUNDS_ZERO_INDEX; //initialized in bounds_initialize
static int    bounds_max_index    = BOUNDS_ZERO_INDEX; //initialized in bounds_initialize

/**
 * @brief Get indexed double precision bin width
 *
 * @return bin width (in bits)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int diwidth() {
  return BIN_WIDTH;
}

/**
 * @brief Get indexed double precision deposit capacity
 *
 * The number of deposits that can be performed before a renorm is necessary. This function applies also to indexed complex double precision.
 *
 * @return deposit capacity
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
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

  bounds[BOUNDS_ZERO_INDEX] = 1.5;
  step = ldexp(1, BIN_WIDTH);

  exp = -1;
  index = BOUNDS_ZERO_INDEX + 1;
  while (exp * BIN_WIDTH  >= DBL_MIN_EXP) {
    bounds[index] = bounds[index - 1] / step;
    index++;
    exp--;
  }
  bounds_max_index = index;
  while (index < BOUNDS_SIZE) {
    bounds[index] = 0.0;
    index++;
  }

  exp = 1;
  index = BOUNDS_ZERO_INDEX - 1;
  while (exp * BIN_WIDTH <= DBL_MAX_EXP) {
    bounds[index] = bounds[index + 1] * step;
    index--;
    exp++;
  }
  bounds_min_index = index;
  while (index >= 0) {
    bounds[index] = bounds[bounds_min_index + 1] * step;
    index--;
  }

  bounds_initialized = 1;
}


/**
 * @internal
 * @brief Get index of manually specified indexed double precision
 *
 * The index of an indexed type is the bin that it corresponds to. Higher indicies correspond to smaller bins.
 *
 * @param repX X's rep vector
 * @return X's index
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int dmindex(const double *repX){
  int index;

  bounds_initialize();

  if(isinf(repX[0])){
    index = bounds_min_index;
  } else if(repX[0] == 0){
    index = bounds_max_index;
  } else {
    frexp(repX[0], &index);
    index--;
    index /= BIN_WIDTH;
    index = BOUNDS_ZERO_INDEX - index;
  }
  return index;
}

/**
 * @brief Get index of double precision
 *
 * The index of a non-indexed type is the smallest index an indexed type would need to have to sum it reproducibly. Higher indicies correspond to smaller bins.
 *
 * @param X scalar X
 * @return X's index
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int dindex(const double X){
  int index;

  bounds_initialize();

  if(isinf(X)){
    index = bounds_min_index;
  }else if(X == 0){
    index = bounds_max_index;
  }else{
    frexp(X, &index);
    index += PREC - BIN_WIDTH;
    if(index < 0){
      index -= BIN_WIDTH - 1; //we want to round towards -infinity
    }
    index /= BIN_WIDTH;
    index = BOUNDS_ZERO_INDEX - 1 - index;
  }
  return index;
}

/**
 * @brief Get double precision bound corresponding to index
 *
 * @param index index
 * @return bound (bin)
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
double dbound(const int index){
  bounds_initialize();

  return bounds[index];
}

/**
 * @internal
 * @brief Set manually specified indexed double precision bounds
 *
 * Set the manually specified indexed double precision X to be empty with the given index
 *
 * @param index index
 * @param repX X's rep vector
 * @param increpX stride within X's rep vector (use every increpX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dmbound(const int fold, const int index, double *repX, const int increpX, double *carX, const int inccarX) {
  int i;

  bounds_initialize();

  for (i = 0; i < fold; i++) {
    repX[i * increpX] = bounds[index + i];
    carX[i * inccarX] = 0.0;
  }
}
