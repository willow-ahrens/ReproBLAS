#include <math.h>

#include <indexed.h>

/**
 * @internal
 * @brief Check if index of manually specified indexed single precision is 0
 *
 * A quick check to determine if the index is 0
 *
 * @param manX X's mantissa vector
 * @return 1 if x has index 0, 0 otherwise.
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
int smindex0(const float *manX){
  int exp;

  frexpf(manX[0], &exp);
  if(exp == FLT_MAX_EXP){
    return 1;
  }
  return 0;
}
