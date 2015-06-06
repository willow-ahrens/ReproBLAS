#include <math.h>

#include <indexed.h>

#include "../../config.h"

static float bins[(FLT_MAX_EXP - FLT_MIN_EXP)/SIWIDTH + MAX_FOLD]; //initialized in bins_initialize
static int   bins_initialized  = 0;                                    //initialized in bins_initialize

static void bins_initialize() {
  int index;

  if (bins_initialized) {
    return;
  }

  bins[0] = ldexpf(0.75, FLT_MAX_EXP);
  for(index = 1; index <= (FLT_MAX_EXP - FLT_MIN_EXP)/SIWIDTH; index++){
    bins[index] = ldexpf(0.75, (FLT_MAX_EXP + FLT_MANT_DIG - SIWIDTH + 1 - index * SIWIDTH));
  }
  for(; index < (FLT_MAX_EXP - FLT_MIN_EXP)/SIWIDTH + MAX_FOLD; index++){
    bins[index] = bins[index - 1];
  }

  bins_initialized = 1;
}

/**
 * @brief Get single precision bin corresponding to index
 *
 * @param X index
 * @return bin (bin)
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
float sbin(const int X){
  bins_initialize();

  return bins[X];
}

/**
 * @internal
 * @brief Set manually specified indexed single precision bins
 *
 * Set the manually specified indexed single precision X to be empty with the given index
 *
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
void smbin(const int fold, const int X, float *manY, const int incmanY, float *carY, const int inccarY) {
  int i;

  bins_initialize();

  for (i = 0; i < fold; i++) {
    manY[i * incmanY] = bins[X + i];
    carY[i * inccarY] = 0.0;
  }
}
