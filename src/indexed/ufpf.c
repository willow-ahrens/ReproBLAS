#include <math.h>

#include "indexed.h"

#include "../common/common.h"

/**
 * @internal
 * @brief unit in the first place
 *
 * This method returns just the implicit 1 in the mantissa of a @c float
 *
 * @param X scalar X
 * @return unit in the first place
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
float ufpf(float X) {
  /*
  int exp;
  if (X == 0.0){
    return 0.0;
  }
  frexpf(X, &exp);
  return ldexpf(0.5, exp);
  */
  int_float x;
  x.f = X;
  x.i &= 0xFF800000ul;
  return x.f;
}
