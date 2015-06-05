#include <indexed.h>

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
