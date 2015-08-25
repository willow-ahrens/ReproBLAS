#include <math.h>

#include <idxd.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Get a reproducible single precision scale
 *
 * For a given X, the smallest Y such that #idxd_sindex(X) == #idxd_sindex(Y)
 *
 * Perhaps the most useful property of this number is that if the bin epsilon of X's bin is normalized, 1.0 <= X * (1.0/Y) < 2^#SIWIDTH
 *
 * @param X single precision number to be scaled
 * @return reproducible scaling factor
 *
 * @author Peter Ahrens
 * @date   19 Jun 2015
 */
float idxd_sscale(const float X){
  return ldexpf(0.5, (FLT_MAX_EXP - SIWIDTH + 1) - idxd_sindex(X) * SIWIDTH);
}
