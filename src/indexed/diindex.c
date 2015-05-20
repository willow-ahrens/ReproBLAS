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
static double bounds[(DBL_MAX_EXP - DBL_MIN_EXP)/BIN_WIDTH + MAX_FOLD + 1]; //initialized in bounds_initialize
static int    bounds_initialized  = 0;                                      //initialized in bounds_initialize

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

/**
 * @internal
 * @brief Get indexed double precision compression factor
 *
 * This factor is used to scale down inputs before deposition
 *
 * @return compression factor
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
double dmcompression() {
  return 1.0/(1 << (PREC - BIN_WIDTH + 1));
}

/**
 * @internal
 * @brief Get indexed double precision expansion factor
 *
 * This factor is used to scale up inputs after deposition
 *
 * @return expansion factor
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
double dmexpansion() {
  return 1.0*(1 << (PREC - BIN_WIDTH + 1));
}

/**
 * @internal
 * @brief Get index of manually specified indexed double precision
 *
 * The index of an indexed type is the bin that it corresponds to. Higher indicies correspond to smaller bins.
 *
 * @param manX X's mantissa vector
 * @return X's index
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   19 May 2015
 */
int dmindex(const double *manX){
  int exp;

  if(manX[0] == 0.0){
    return (DBL_MAX_EXP - DBL_MIN_EXP)/BIN_WIDTH + 1;
  }else{
    frexp(manX[0], &exp);
    return (DBL_MAX_EXP - exp)/BIN_WIDTH;
  }
}

/**
 * @brief Get index of double precision
 *
 * The index of a non-indexed type is the smallest index an indexed type would need to have to sum it reproducibly. Higher indicies correspond to smaller bins.
 *
 * @param X scalar X
 * @return X's index
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   19 May 2015
 */
int dindex(const double X){
  int exp;

  if(X == 0.0){
    return (DBL_MAX_EXP - DBL_MIN_EXP)/BIN_WIDTH + 1;
  }else{
    frexp(X, &exp);
    return (DBL_MAX_EXP - exp)/BIN_WIDTH;
  }
}

static void bounds_initialize() {
  int index;

  if (bounds_initialized) {
    return;
  }

  for(index = 0; index <= (DBL_MAX_EXP - DBL_MIN_EXP)/BIN_WIDTH; index++){
    bounds[index] = ldexp(0.75, (DBL_MAX_EXP - index * BIN_WIDTH));
  }
  for(; index < (DBL_MAX_EXP - DBL_MIN_EXP)/BIN_WIDTH + MAX_FOLD + 1; index++){
    bounds[index] = 0;
  }

  bounds_initialized = 1;
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
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dmbound(const int fold, const int index, double *manX, const int incmanX, double *carX, const int inccarX) {
  int i;

  bounds_initialize();

  for (i = 0; i < fold; i++) {
    manX[i * incmanX] = bounds[index + i];
    carX[i * inccarX] = 0.0;
  }
}
