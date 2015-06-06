#include <math.h>

#include <indexed.h>

#include "../../config.h"

static double bins[(DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH + MAX_FOLD]; //initialized in bins_initialize
static int    bins_initialized  = 0;                                    //initialized in bins_initialize

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
