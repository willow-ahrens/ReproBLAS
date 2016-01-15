#include <math.h>

#include <idxd.h>

#include "../../config.h"

static float bins[idxd_SIMAXINDEX + idxd_SIMAXFOLD];
static int bins_initialized = 0;

/**
 * @internal
 * @brief Get indexed single precision reference bins
 *
 * returns a pointer to the bins corresponding to the given index
 *
 * @param X index
 * @return pointer to constant single precision bins of index X
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   19 Jun 2015
 */
const float *idxd_smbins(const int X) {
  int index;

  if (!bins_initialized){
    bins[0] = ldexpf(0.75, FLT_MAX_EXP);
    for(index = 1; index <= idxd_SIMAXINDEX; index++){
      bins[index] = ldexpf(0.75, (FLT_MAX_EXP + FLT_MANT_DIG - SIWIDTH + 1 - index * SIWIDTH));
    }
    for(; index < idxd_SIMAXINDEX + idxd_SIMAXFOLD; index++){
      bins[index] = bins[index - 1];
    }

    bins_initialized = 1;
  }

  return (const float*)bins + X;
}
