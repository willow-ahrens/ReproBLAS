#include <indexed.h>

/**
 * @internal
 * @brief Negate manually specified indexed single precision (X = -X)
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
void smnegate(const int fold, float* manX, const int incmanX, float* carX, const int inccarX) {
  int i;
  for (i = 0; i < fold; i++) {
    manX[i * incmanX] = (3 * ufp(manX[i * incmanX])) - manX[i * incmanX];
    carX[i * inccarX] = -carX[i * inccarX];
  }
}
