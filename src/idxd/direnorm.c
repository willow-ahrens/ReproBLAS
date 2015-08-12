#include <idxd.h>

/**
 * @brief Renormalize indexed double precision
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
void direnorm(const int fold, double_indexed *X) {
  dmrenorm(fold, X, 1, X + fold, 1);
}
