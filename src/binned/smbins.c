#include <math.h>

#include <binned.h>

#include "../../config.h"

static float bins[binned_SBMAXINDEX + binned_SBMAXFOLD];
static int bins_initialized = 0;

/**
 * @internal
 * @brief Get binned single precision reference bins
 *
 * returns a pointer to the bins corresponding to the given index
 *
 * @param X index
 * @return pointer to constant single precision bins of index X
 *
 * @author Willow Ahrens
 * @author Hong Diep Nguyen
 * @date   19 Jun 2015
 */
const float *binned_smbins(const int X) {
  int index;

  if (!bins_initialized){
    bins[0] = ldexpf(0.75, FLT_MAX_EXP);
    for(index = 1; index <= binned_SBMAXINDEX; index++){
      bins[index] = ldexpf(0.75, (FLT_MAX_EXP + FLT_MANT_DIG - SBWIDTH + 1 - index * SBWIDTH));
    }
    for(; index < binned_SBMAXINDEX + binned_SBMAXFOLD; index++){
      bins[index] = bins[index - 1];
    }

    bins_initialized = 1;
  }

  return (const float*)bins + X;
}
