#include <math.h>

#include <indexed.h>

#include "../../config.h"
#include "../common/common.h"

#define DIWIDTH        40
static double bins[(DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH + MAX_FOLD]; //initialized in bins_initialize
static int    bins_initialized  = 0;                                    //initialized in bins_initialize

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
  return X * ldexp(0.5, (1 - fold)*(DIWIDTH - 1) + 1) * N;
}

/**
 * @internal
 * @brief Get a reproducible double precision scale
 *
 * The scaling factor Y returned for given X is the smallest value that will fit in X's bin (The smallest representable value with the same index as X)
 *
 * Perhaps the most useful property of this number is that 1.0 <= X * (1.0/Y) < 2^#DIWIDTH
 *
 * @param X double precision number to be scaled
 * @return reproducible scaling factor (if X == 0.0, returns smallest valid scale)
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
double dscale(const double X){
  return ldexp(0.5, (DBL_MAX_EXP - DIWIDTH + 1 - MIN(dindex(X), (DBL_MAX_EXP - DBL_MIN_EXP - DIWIDTH)/DIWIDTH) * DIWIDTH));
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
    return (DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH + MAX_FOLD;
  }else{
    frexp(manX[0], &exp);
    if(exp == DBL_MAX_EXP){
      return 0;
    }
    return (DBL_MAX_EXP + DBL_MANT_DIG - DIWIDTH + 1 - exp)/DIWIDTH;
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
    return (DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH;
  }else{
    frexp(X, &exp);
    return (DBL_MAX_EXP - exp)/DIWIDTH;
  }
}

static void bins_initialize() {
  int index;

  if (bins_initialized) {
    return;
  }

  bins[0] = ldexp(0.75, DBL_MAX_EXP);
  for(index = 1; index <= (DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH; index++){
    bins[index] = ldexp(0.75, (DBL_MAX_EXP + DBL_MANT_DIG - DIWIDTH + 1 - index * DIWIDTH));
  }
  for(; index < (DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH + MAX_FOLD; index++){
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
