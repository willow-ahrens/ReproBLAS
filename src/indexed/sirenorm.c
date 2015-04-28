/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

/**
 * @internal
 * @brief Renormalize manually specified indexed single precision
 *
 * Renormalization keeps the rep vector within the necessary bounds by shifting over to the carry vector
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
void smrenorm(const int fold, float* repX, int increpX, float* carX, int inccarX) {
  int i;
  float M;
  float repX0;
  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    repX0 = repX[0];
    if (repX0 == 0.0)
      continue;

    M = ufpf(repX0);
    if (repX0 >= (M * 1.75)) {
      repX[0] -= M * 0.25;
      carX[0] += 1;
    }
    else if (repX0 < (M * 1.25)) {
      repX[0] += M * 0.5;
      carX[0] -= 2;
    }
    else if (repX0 < (M * 1.5)) {
      repX[0] += M * 0.25;
      carX[0] -= 1;
    }
  }
}

/**
 * @brief Renormalize indexed single precision
 *
 * Renormalization keeps the rep vector within the necessary bounds by shifting over to the carry vector
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void sirenorm(const int fold, float_indexed *X) {
  smrenorm(fold, X, 1, X + fold, 1);
}

/**
 * @internal
 * @brief Renormalize manually specified indexed complex single precision
 *
 * Renormalization keeps the rep vector within the necessary bounds by shifting over to the carry vector
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
void cmrenorm(const int fold, float* repX, int increpX, float* carX, int inccarX) {
  smrenorm(fold, repX, 2 * increpX, carX, 2 * inccarX);
  smrenorm(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX);
}

/**
 * @brief Renormalize indexed complex single precision
 *
 * Renormalization keeps the rep vector within the necessary bounds by shifting over to the carry vector
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cirenorm(const int fold, float_complex_indexed *X) {
  cmrenorm(fold, X, 1, X + 2 * fold, 1);
}
