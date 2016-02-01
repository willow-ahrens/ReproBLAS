#include <idxd.h>

/**
 * @brief Print indexed double precision
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void idxd_diprint(const int fold, const double_indexed *X) {
  idxd_dmprint(fold, X, 1, X + fold, 1);
}
