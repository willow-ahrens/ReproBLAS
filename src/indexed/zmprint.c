#include <math.h>
#include <stdio.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Print manually specified indexed complex double precision
 *
 * @param fold the fold of the indexed types
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zmprint(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX){
  int i;
  double M;
  for (i = 0; i < fold; i++, priX += incpriX, carX += inccarX) {
    M = UFP(priX[0]);
    printf("(2^%d: %g #%g =%g", (int)log2(M) + 1, priX[0] - 1.5*M, carX[0], ((carX[0] - 6) * 0.25 * M + priX[0]));
    M = UFP(priX[1]);
    printf("| 2^%d: %g #%g =%g)\n", (int)log2(M) + 1, priX[1] - 1.5*M, carX[1], ((carX[1] - 6) * 0.25 * M + priX[1]));
  }
}
