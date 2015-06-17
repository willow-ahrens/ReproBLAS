#include <math.h>
#include <stdio.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Print manually specified indexed complex double precision
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
void zmprint(const int fold, const double *manX, const int incmanX, const double *carX, const int inccarX){
  int i;
  double M;
  for (i = 0; i < fold; i++, manX += incmanX, carX += inccarX) {
    M = UFP(manX[0]);
    printf("(2^%d: %g #%g =%g", (int)log2(M) + 1, manX[0] - 1.5*M, carX[0], ((carX[0] - 6) * 0.25 * M + manX[0]));
    M = UFP(manX[1]);
    printf("| 2^%d: %g #%g =%g)\n", (int)log2(M) + 1, manX[1] - 1.5*M, carX[1], ((carX[1] - 6) * 0.25 * M + manX[1]));
  }
}
