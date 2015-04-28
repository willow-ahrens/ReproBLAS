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
 * @brief Renormalize manually specified indexed double precision
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
void dmrenorm(const int fold, double* repX, const int increpX, double* carX, const int inccarX) {
  int i;
  double M;
  double repX0;
  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    repX0 = repX[0];
    if (repX0 == 0.0)
      continue;

    M = ufp(repX0);
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
 * @brief Renormalize indexed double precision
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
void direnorm(const int fold, double_indexed *X) {
  dmrenorm(fold, X, 1, X + fold, 1);
}

/**
 * @internal
 * @brief Renormalize manually specified indexed complex double precision
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
void zmrenorm(const int fold, double* repX, const int increpX, double* carX, const int inccarX) {
  dmrenorm(fold, repX, 2 * increpX, carX, 2 * inccarX);
  dmrenorm(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX);
}

/**
 * @brief Renormalize indexed complex double precision
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
void zirenorm(const int fold, double_complex_indexed *X) {
  zmrenorm(fold, X, 1, X + 2 * fold, 1);
}
