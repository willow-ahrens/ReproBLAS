#include "indexed.h"

/**
 * @internal
 * @brief Negate manually specified indexed double precision (X = -X)
 *
 * Performs the operation X = -X
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
void dmnegate(const int fold, double* manX, const int incmanX, double* carX, const int inccarX) {
  int i;
  for (i = 0; i < fold; i++) {
    manX[i * incmanX] = (3 * ufp(manX[i * incmanX])) - manX[i * incmanX];
    carX[i * inccarX] = -carX[i * inccarX];
  }
}

/**
 * @brief Negate indexed double precision (X = -X)
 *
 * Performs the operation X = -X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dinegate(const int fold, double_indexed* X){
  dmnegate(fold, X, 1, X + fold, 1);
}

/**
 * @internal
 * @brief Negate manually specified indexed complex double precision (X = -X)
 *
 * Performs the operation X = -X
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
void zmnegate(const int fold, double* manX, const int incmanX, double* carX, const int inccarX) {
  dmnegate(fold, manX, 2 * incmanX, carX, 2 * inccarX);
  dmnegate(fold, manX + 1, 2 * incmanX, carX + 1, 2 * inccarX);
}

/**
 * @brief Negate indexed complex double precision (X = -X)
 *
 * Performs the operation X = -X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zinegate(const int fold, double_complex_indexed* X){
  zmnegate(fold, X, 1, X + 2 * fold, 1);
}
