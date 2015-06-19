#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Get a reproducible single precision scale
 *
 * If bin epsilon is the smallest value that will fit in a given bin, #sscale(X) returns Y, the smallest normalized bin epsilon such that #sindex(Y) >= #sindex(X)
 *
 * Perhaps the most useful property of this number is that if the bin epsilon of X's bin is normalized, 1.0 <= X * (1.0/Y) < 2^#SIWIDTH
 *
 * @param X single precision number to be scaled
 * @return reproducible scaling factor
 *
 * @author Peter Ahrens
 * @date   19 Jun 2015
 */
float sscale(const float X){
  return ldexpf(0.5, (FLT_MAX_EXP - SIWIDTH + 1 - MIN(sindex(X), (FLT_MAX_EXP - FLT_MIN_EXP - SIWIDTH)/SIWIDTH) * SIWIDTH));
}
