#include <idxd.h>

/**
 * @brief Print indexed single precision
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void siprint(const int fold, const float_indexed *X) {
  smprint(fold, X, 1, X + fold, 1);
}
