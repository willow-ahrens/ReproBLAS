/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <math.h>
#include <float.h>

#include "indexed.h"
#include "../../config.h"
#include "../common/common.h"

#define DBL_BIN_DIG        40
static double bins[(DBL_MAX_EXP - DBL_MIN_EXP)/DBL_BIN_DIG + MAX_FOLD]; //initialized in bins_initialize
static int    bins_initialized  = 0;                                    //initialized in bins_initialize

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
  return DBL_BIN_DIG;
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
  return 1 << (DBL_MANT_DIG - DBL_BIN_DIG - 2);
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
  return 1.0/(1 << (DBL_MANT_DIG - DBL_BIN_DIG + 1));
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
  return 1.0*(1 << (DBL_MANT_DIG - DBL_BIN_DIG + 1));
}

/**
 * @brief Get indexed double precision summation error bound
 *
 * This is a bound on the absolute error of a summation using indexed types
 *
 * @param fold the fold of the indexed types
 * @param N the number of double precision floating point summands
 * @param X the maximum absolute value of the summands
 * @return error bound
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   21 May 2015
 */
double dibound(const int fold, const int N, const double X) {
  return X * ldexp(0.5, (1 - fold)*(DBL_BIN_DIG - 1) + 1) * N;
}

/**
 * @internal
 * @brief Get a reproducible double precision scale
 *
 * The scaling factor Y returned for given X is the smallest value that will fit in X's bin (The smallest representable value with the same index as X)
 *
 * Perhaps the most useful property of this number is that 1.0 <= X * (1.0/Y) < 2^#diwidth()
 *
 * @param X double precision number to be scaled
 * @return reproducible scaling factor (if X == 0.0, returns smallest valid scale)
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
double dscale(const double X){
  return ldexp(0.5, MAX((DBL_MAX_EXP - DBL_BIN_DIG + 1 - dindex(X) * DBL_BIN_DIG), DBL_MIN_EXP));
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
    return (DBL_MAX_EXP - DBL_MIN_EXP)/DBL_BIN_DIG + MAX_FOLD;
  }else{
    frexp(manX[0], &exp);
    if(exp == DBL_MAX_EXP){
      return 0;
    }
    return (DBL_MAX_EXP + DBL_MANT_DIG - DBL_BIN_DIG + 1 - exp)/DBL_BIN_DIG;
  }
}

/**
 * @internal
 * @brief Check if index of manually specified indexed double precision is 0
 *
 * A quick check to determine if the index is 0
 *
 * @param manX X's mantissa vector
 * @return 1 if x has index 0, 0 otherwise.
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
int dmindex0(const double *manX){
  int exp;

  frexp(manX[0], &exp);
  if(exp == DBL_MAX_EXP){
    return 1;
  }
  return 0;
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
    return (DBL_MAX_EXP - DBL_MIN_EXP)/DBL_BIN_DIG;
  }else{
    frexp(X, &exp);
    return (DBL_MAX_EXP - exp)/DBL_BIN_DIG;
  }
}

static void bins_initialize() {
  int index;

  if (bins_initialized) {
    return;
  }

  bins[0] = ldexp(0.75, DBL_MAX_EXP);
  for(index = 1; index <= (DBL_MAX_EXP - DBL_MIN_EXP)/DBL_BIN_DIG; index++){
    bins[index] = ldexp(0.75, (DBL_MAX_EXP + DBL_MANT_DIG - DBL_BIN_DIG + 1 - index * DBL_BIN_DIG));
  }
  for(; index < (DBL_MAX_EXP - DBL_MIN_EXP)/DBL_BIN_DIG + MAX_FOLD; index++){
    bins[index] = bins[index - 1];
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
