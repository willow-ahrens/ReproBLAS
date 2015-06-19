#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Get a reproducible double precision scale
 *
 * If bin epsilon is the smallest value that will fit in a given bin, #dscale(X) returns Y, the smallest normalized bin epsilon such that #dindex(Y) >= #dindex(X)
 *
 * Perhaps the most useful property of this number is that if the bin epsilon of X's bin is normalized, 1.0 <= X * (1.0/Y) < 2^#DIWIDTH
 *
 * @param X double precision number to be scaled
 * @return reproducible scaling factor
 *
 * @author Peter Ahrens
 * @date   19 Jun 2015
 */
double dscale(const double X){
  return ldexp(0.5, (DBL_MAX_EXP - DIWIDTH + 1 - MIN(dindex(X), (DBL_MAX_EXP - DBL_MIN_EXP - DIWIDTH)/DIWIDTH) * DIWIDTH));
}
