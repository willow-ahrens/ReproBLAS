#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Get a reproducible double precision scale
 *
 * For a given X, the smallest Y such that #dindex(X) == #dindex(Y)
 *
 * Perhaps the most useful property of this number is that, 1.0 <= X/Y < 2^#DIWIDTH
 *
 * @param X double precision number to be scaled
 * @return reproducible scaling factor
 *
 * @author Peter Ahrens
 * @date   19 Jun 2015
 */
double dscale(const double X){
  return ldexp(0.5, (DBL_MAX_EXP - DIWIDTH + 1) - dindex(X) * DIWIDTH);
}
