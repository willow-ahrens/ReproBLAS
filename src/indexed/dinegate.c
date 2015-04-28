#include "indexed.h"

/**
 * @internal
 * @brief Negate manually specified indexed double precision (X = -X)
 *
 * Performs the operation X = -X
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
void dmnegate(const int fold, double* repX, const int increpX, double* carX, const int inccarX) {
  int i;
  for (i = 0; i < fold; i++) {
    repX[i * increpX] = (3 * ufp(repX[i * increpX])) - repX[i * increpX];
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
 * @param repX X's rep vector
 * @param increpX stride within X's rep vector (use every increpX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zmnegate(const int fold, double* repX, const int increpX, double* carX, const int inccarX) {
  dmnegate(fold, repX, 2 * increpX, carX, 2 * inccarX);
  dmnegate(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX);
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
