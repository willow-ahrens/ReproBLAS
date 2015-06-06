#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Get a reproducible single precision scale
 *
 * The scaling factor Y returned for given X is the smallest value that will fit in X's bin (The smallest representable value with the same index as X)
 *
 * Perhaps the most useful property of this number is that 1.0 <= X * (1.0/Y) < 2^#SIWIDTH
 *
 * @param X single precision number to be scaled
 * @return reproducible scaling factor (if X == 0.0, returns smallest valid scale)
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
float sscale(const float X){
  return ldexpf(0.5, (FLT_MAX_EXP - SIWIDTH + 1 - MIN(sindex(X), (FLT_MAX_EXP - FLT_MIN_EXP - SIWIDTH)/SIWIDTH) * SIWIDTH));
}
