#include <math.h>

#include "indexed.h"

double ufp(double x) {
  int exp;
  if (x == 0.0) {
    return 0.0;
  }
  frexp(x, &exp);
  return ldexp(0.5, exp);
}
