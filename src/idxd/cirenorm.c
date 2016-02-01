#include <idxd.h>

/**
 * @brief Renormalize indexed complex single precision
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
void idxd_cirenorm(const int fold, float_complex_indexed *X) {
  idxd_cmrenorm(fold, X, 1, X + 2 * fold, 1);
}
