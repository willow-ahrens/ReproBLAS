#include <indexed.h>

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
  const double *bins;
  if(manX[0] != 0.0){
    bins = dmbins(dmindex(manX));
    for (i = 0; i < fold; i++) {
      manX[i * incmanX] = bins[i] - (manX[i * incmanX] - bins[i]);
      carX[i * inccarX] = -carX[i * inccarX];
    }
  }
}
