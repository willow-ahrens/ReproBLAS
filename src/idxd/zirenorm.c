#include <idxd.h>

/**
 * @brief Renormalize indexed complex double precision
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
void idxd_zirenorm(const int fold, double_complex_indexed *X) {
  idxd_zmrenorm(fold, X, 1, X + 2 * fold, 1);
}
