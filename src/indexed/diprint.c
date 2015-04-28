#include <stdio.h>

#include "indexed.h"

/**
 * @internal
 * @brief Print manually specified indexed double precision
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
void dmprint(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX) {
  int i;
  double M;
  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    M = ufp(repX[0]);
    printf("(2^%d: %g #%g =%g)\n", (int)log2(M), repX[0] - 1.5*M, carX[0], (carX[0] - 6) * 0.25 * M + repX[0]);
  }
}

/**
 * @brief Print indexed double precision
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void diprint(const int fold, const double_indexed *X) {
  dmprint(fold, X, 1, X + fold, 1);
}

/**
 * @internal
 * @brief Print manually specified indexed complex double precision
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
void zmprint(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX){
  int i;
  double M;
  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    M = ufp(repX[0]);
    printf("(2^%d: %g #%g =%g", (int)log2(M), repX[0] - 1.5*M, carX[0], (carX[0] - 6) * 0.25 * M + repX[0]);
    M = ufp(repX[1]);
    printf("| 2^%d: %g #%g =%g)\n", (int)log2(M), repX[1] - 1.5*M, carX[1], (carX[1] - 6) * 0.25 * M + repX[1]);
  }
}

/**
 * @brief Print indexed complex double precision
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void ziprint(const int fold, const double_complex_indexed *X){
  zmprint(fold, X, 2, X + 2 * fold, 2);
}
