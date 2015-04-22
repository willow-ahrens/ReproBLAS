#include <math.h>

#include "indexed.h"

float ufpf(float x) {
  int exp;
  if (x == 0.0){
    return 0.0;
  }
  frexpf(x, &exp);
  return ldexpf(0.5, exp);
}
