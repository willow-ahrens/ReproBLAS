#include "indexed.h"

/**
 * @internal
 * @brief Renormalize manually specified indexed double precision
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
void dmrenorm(const int fold, double* manX, const int incmanX, double* carX, const int inccarX) {
  int i;
  double M;
  double manX0;
  for (i = 0; i < fold; i++, manX += incmanX, carX += inccarX) {
    manX0 = manX[0];
    if (manX0 == 0.0)
      continue;

    M = ufp(manX0);
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
 * @brief Renormalize indexed double precision
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
void direnorm(const int fold, double_indexed *X) {
  dmrenorm(fold, X, 1, X + fold, 1);
}

/**
 * @internal
 * @brief Renormalize manually specified indexed complex double precision
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
void zmrenorm(const int fold, double* manX, const int incmanX, double* carX, const int inccarX) {
  dmrenorm(fold, manX, 2 * incmanX, carX, 2 * inccarX);
  dmrenorm(fold, manX + 1, 2 * incmanX, carX + 1, 2 * inccarX);
}

/**
 * @brief Renormalize indexed complex double precision
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
void zirenorm(const int fold, double_complex_indexed *X) {
  zmrenorm(fold, X, 1, X + 2 * fold, 1);
}
