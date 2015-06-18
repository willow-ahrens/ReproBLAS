#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Renormalize manually specified indexed single precision
 *
 * Renormalization keeps the mantissa vector within the necessary bins by shifting over to the carry vector
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void smrenorm(const int fold, float* manX, const int incmanX, float* carX, const int inccarX) {
  int i;
  float M;
  float manX0 = manX[0];

  if(manX0 == 0.0 || ISNANINFF(manX0)){
    return;
  }

  for (i = 0; i < fold; i++, manX += incmanX, carX += inccarX) {
    manX0 = manX[0];

    M = UFPF(manX0);

    if (manX0 >= (M * 1.75)) {
      manX[0] -= M * 0.25;
      carX[0] += 1;
    }
    else if (manX0 < (M * 1.25)) {
      manX[0] += M * 0.5;
      carX[0] -= 2;
    }
    else if (manX0 < (M * 1.5)) {
      manX[0] += M * 0.25;
      carX[0] -= 1;
    }
  }
}
