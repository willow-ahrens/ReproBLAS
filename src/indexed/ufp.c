#include <math.h>

#include "indexed.h"

#include "../common/common.h"

/**
 * @internal
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
double ufp(double X) {
  /*
  int exp;
  if (X == 0.0) {
    return 0.0;
  }
  frexp(X, &exp);
  return ldexp(0.5, exp);
  */
  //faster
  long_double x;
  x.d = X;
  x.l &= 0xFFF0000000000000ull;
  return x.d;
}
