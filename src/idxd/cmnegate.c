#include <idxd.h>

/**
 * @internal
 * @brief Negate manually specified indexed complex single precision (X = -X)
 *
 * Performs the operation X = -X
 *
 * @param fold the fold of the indexed types
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmnegate(const int fold, float* priX, const int incpriX, float* carX, const int inccarX) {
  smnegate(fold, priX, 2 * incpriX, carX, 2 * inccarX);
  smnegate(fold, priX + 1, 2 * incpriX, carX + 1, 2 * inccarX);
}
