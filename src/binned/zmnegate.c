#include <binned.h>

/**
 * @internal
 * @brief Negate manually specified binned complex double precision (X = -X)
 *
 * Performs the operation X = -X
 *
 * @param fold the fold of the binned types
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_zmnegate(const int fold, double* priX, const int incpriX, double* carX, const int inccarX) {
  binned_dmnegate(fold, priX, 2 * incpriX, carX, 2 * inccarX);
  binned_dmnegate(fold, priX + 1, 2 * incpriX, carX + 1, 2 * inccarX);
}
