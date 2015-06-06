#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Get a reproducible double precision scale
 *
 * The scaling factor Y returned for given X is the smallest value that will fit in X's bin (The smallest representable value with the same index as X)
 *
 * Perhaps the most useful property of this number is that 1.0 <= X * (1.0/Y) < 2^#DIWIDTH
 *
 * @param X double precision number to be scaled
 * @return reproducible scaling factor (if X == 0.0, returns smallest valid scale)
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
double dscale(const double X){
  return ldexp(0.5, (DBL_MAX_EXP - DIWIDTH + 1 - MIN(dindex(X), (DBL_MAX_EXP - DBL_MIN_EXP - DIWIDTH)/DIWIDTH) * DIWIDTH));
}
