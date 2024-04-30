#include <math.h>
#include <stdio.h>

#include <binned.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Print manually specified binned complex single precision
 *
 * @param fold the fold of the binned types
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_cmprint(const int fold, const float* priX, const int incpriX, const float* carX, const int inccarX) {
  int i;
  float M;
  for (i = 0; i < fold; i++, priX += incpriX, carX += inccarX) {
    M = UFPF(priX[0]);
    printf("(2^%d: %g #%g =%g", (int)log2f(M) + 1, priX[0] - 1.5*M, carX[0], ((carX[0] - 6) * 0.25 * M + priX[0]));
    M = UFPF(priX[1]);
    printf("| 2^%d: %g #%g =%g)\n", (int)log2f(M) + 1, priX[1] - 1.5*M, carX[1], ((carX[1] - 6) * 0.25 * M + priX[1]));
  }
}
