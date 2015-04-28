#include <stdio.h>

#include "indexed.h"

/**
 * @internal
 * @brief Print manually specified indexed single precision
 *
 * @param fold the fold of the indexed types
 * @param repX X's rep vector
 * @param increpX stride within X's rep vector (use every increpX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void smprint(const int fold, const float *repX, const int increpX, const float *carX, const int inccarX) {
  int i;
  float M;
  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    M = ufpf(repX[0]);
    printf("(2^%d: %g #%g =%g)\n", (int)log2f(M), repX[0] - 1.5*M, carX[0], (carX[0] - 6) * 0.25 * M + repX[0]);
  }
}

/**
 * @brief Print indexed single precision
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void siprint(const int fold, const float_indexed *X) {
  smprint(fold, X, 1, X + fold, 1);
}

/**
 * @internal
 * @brief Print manually specified indexed complex single precision
 *
 * @param fold the fold of the indexed types
 * @param repX X's rep vector
 * @param increpX stride within X's rep vector (use every increpX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmprint(const int fold, const float* repX, const int increpX, const float* carX, const int inccarX) {
  int i;
  float M;
  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    M = ufpf(repX[0]);
    printf("(2^%d: %g #%g =%g", (int)log2f(M), repX[0] - 1.5*M, carX[0], (carX[0] - 6) * 0.25 * M + repX[0]);
    M = ufpf(repX[1]);
    printf("| 2^%d: %g #%g =%g)\n", (int)log2f(M), repX[1] - 1.5*M, carX[1], (carX[1] - 6) * 0.25 * M + repX[1]);
  }
}

/**
 * @brief Print indexed complex single precision
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void ciprint(const int fold, const float_complex_indexed *X) {
  cmprint(fold, X, 2, X + fold * 2, 2);
}
