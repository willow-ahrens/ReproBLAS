#include <math.h>

#include "indexed.h"

/**
 * @brief unit in the first place
 *
 * This method returns just the implicit 1 in the mantissa of a @c double
 *
 * @param X scalar X
 * @return unit in the first place
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
double ufp(double x) {
  int exp;
  if (x == 0.0) {
    return 0.0;
  }
  frexp(x, &exp);
  return ldexp(0.5, exp);
}
