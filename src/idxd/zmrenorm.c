#include <idxd.h>

/**
 * @internal
 * @brief Renormalize manually specified indexed complex double precision
 *
 * Renormalization keeps the primary vector within the necessary bins by shifting over to the carry vector
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
void zmrenorm(const int fold, double* priX, const int incpriX, double* carX, const int inccarX) {
  dmrenorm(fold, priX, 2 * incpriX, carX, 2 * inccarX);
  dmrenorm(fold, priX + 1, 2 * incpriX, carX + 1, 2 * inccarX);
}
