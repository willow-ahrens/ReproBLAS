#include <idxd.h>

/**
 * @brief Renormalize indexed single precision
 *
 * Renormalization keeps the primary vector within the necessary bins by shifting over to the carry vector
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void idxd_sirenorm(const int fold, float_indexed *X) {
  idxd_smrenorm(fold, X, 1, X + fold, 1);
}
