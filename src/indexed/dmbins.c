#include <math.h>

#include <indexed.h>

#include "../../config.h"

static double bins[(DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH + MAX_FOLD];
static int bins_initialized = 0;

/**
 * @internal
 * @brief Get indexed double precision reference bins
 *
 * returns a pointer to the bins corresponding to the given index
 *
 * @param X index
 * @return pointer to constant double precision bins of index X
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   5 Jun 2015
 */
const double *dmbins(const int X) {
  int index;

  if (!bins_initialized) {
    bins[0] = ldexp(0.75, DBL_MAX_EXP);
    for(index = 1; index <= (DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH; index++){
      bins[index] = ldexp(0.75, (DBL_MAX_EXP + DBL_MANT_DIG - DIWIDTH + 1 - index * DIWIDTH));
    }
    for(; index < (DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH + MAX_FOLD; index++){
      bins[index] = bins[index - 1];
    }

    bins_initialized = 1;
  }

  return (const double*)bins + X;
}
