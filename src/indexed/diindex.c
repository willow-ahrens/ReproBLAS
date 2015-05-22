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
static double bins[(DBL_MAX_EXP - DBL_MIN_EXP)/BIN_WIDTH + MAX_FOLD + 1]; //initialized in bins_initialize
static int    bins_initialized  = 0;                                      //initialized in bins_initialize

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
 * @brief Get indexed double precision summation error bin
 *
 * This is a bin on the absolute error of a summation using indexed types
 *
 * @param fold the fold of the indexed types
 * @param N the number of double precision floating point summands
 * @param X the maximum absolute value of the summands
 * @return error bin
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   21 May 2015
 */
double dibound(const int fold, const int N, const double X) {
  return X * ldexp(0.5, (1 - fold)*(BIN_WIDTH - 1) + 1) * N;
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

static void bins_initialize() {
  int index;

  if (bins_initialized) {
    return;
  }

  for(index = 0; index <= (DBL_MAX_EXP - DBL_MIN_EXP)/BIN_WIDTH; index++){
    bins[index] = ldexp(0.75, (DBL_MAX_EXP - index * BIN_WIDTH));
  }
  for(; index < (DBL_MAX_EXP - DBL_MIN_EXP)/BIN_WIDTH + MAX_FOLD + 1; index++){
    bins[index] = 0;
  }

  bins_initialized = 1;
}


/**
 * @brief Get double precision bin corresponding to index
 *
 * @param X index
 * @return bin
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
double dbin(const int X){
  bins_initialize();

  return bins[X];
}

/**
 * @internal
 * @brief Set manually specified indexed double precision bins
 *
 * Set the manually specified indexed double precision X to be empty with the given index
 *
 * @param fold the fold of the indexed type
 * @param X index
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dmbin(const int fold, const int X, double *manY, const int incmanY, double *carY, const int inccarY) {
  int i;

  bins_initialize();

  for (i = 0; i < fold; i++) {
    manY[i * incmanY] = bins[X + i];
    carY[i * inccarY] = 0.0;
  }
}
