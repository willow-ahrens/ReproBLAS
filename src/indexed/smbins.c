#include <math.h>

#include <indexed.h>

#include "../../config.h"

static float bins[(FLT_MAX_EXP - FLT_MIN_EXP)/SIWIDTH + MAX_FOLD];
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
 * @date   5 Jun 2015
 */
const float *smbins(const int X) {
  int index;

  if (!bins_initialized){
    bins[0] = ldexpf(0.75, FLT_MAX_EXP);
    for(index = 1; index <= (FLT_MAX_EXP - FLT_MIN_EXP)/SIWIDTH; index++){
      bins[index] = ldexpf(0.75, (FLT_MAX_EXP + FLT_MANT_DIG - SIWIDTH + 1 - index * SIWIDTH));
    }
    for(; index < (FLT_MAX_EXP - FLT_MIN_EXP)/SIWIDTH + MAX_FOLD; index++){
      bins[index] = bins[index - 1];
    }

    bins_initialized = 1;
  }

  return (const float*)bins + X;
}
