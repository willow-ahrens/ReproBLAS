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
  float manX0;
  for (i = 0; i < fold; i++, manX += incmanX, carX += inccarX) {
    manX0 = manX[0];
    if (manX0 == 0.0)
      continue;

    M = ufpf(manX0);
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

/**
 * @brief Renormalize indexed single precision
 *
 * Renormalization keeps the mantissa vector within the necessary bins by shifting over to the carry vector
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
void cmrenorm(const int fold, float* manX, const int incmanX, float* carX, const int inccarX) {
  smrenorm(fold, manX, 2 * incmanX, carX, 2 * inccarX);
  smrenorm(fold, manX + 1, 2 * incmanX, carX + 1, 2 * inccarX);
}

/**
 * @brief Renormalize indexed complex single precision
 *
 * Renormalization keeps the mantissa vector within the necessary bins by shifting over to the carry vector
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
