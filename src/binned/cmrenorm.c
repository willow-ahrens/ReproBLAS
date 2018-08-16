#include <binned.h>

/**
 * @internal
 * @brief Renormalize manually specified binned complex single precision
 *
 * Renormalization keeps the primary vector within the necessary bins by shifting over to the carry vector
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
void binned_cmrenorm(const int fold, float* priX, const int incpriX, float* carX, const int inccarX) {
  binned_smrenorm(fold, priX, 2 * incpriX, carX, 2 * inccarX);
  binned_smrenorm(fold, priX + 1, 2 * incpriX, carX + 1, 2 * inccarX);
}
