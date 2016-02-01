#include <math.h>

#include <idxd.h>

#include "../../config.h"

static double bins[idxd_DIMAXINDEX + idxd_DIMAXFOLD];
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
 * @date   19 Jun 2015
 */
const double *idxd_dmbins(const int X) {
  int index;

  if (!bins_initialized) {
    bins[0] = 2.0 * ldexp(0.75, DBL_MAX_EXP - 1);
    for(index = 1; index <= idxd_DIMAXINDEX; index++){
      bins[index] = ldexp(0.75, (DBL_MAX_EXP + DBL_MANT_DIG - DIWIDTH + 1 - index * DIWIDTH));
    }
    for(; index < idxd_DIMAXINDEX + idxd_DIMAXFOLD; index++){
      bins[index] = bins[index - 1];
    }

    bins_initialized = 1;
  }

  return (const double*)bins + X;
}
