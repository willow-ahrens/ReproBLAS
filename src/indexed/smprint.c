#include <math.h>
#include <stdio.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Print manually specified indexed single precision
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
void smprint(const int fold, const float *manX, const int incmanX, const float *carX, const int inccarX) {
  int i;
  float M;
  for (i = 0; i < fold; i++, manX += incmanX, carX += inccarX) {
    M = UFPF(manX[0]);
    printf("(2^%d: %g #%g =%g)\n", (int)log2f(M) + 1, manX[0] - 1.5*M, carX[0], ((carX[0] - 6) * 0.25 * M + manX[0]));
  }
}
