#include <indexed.h>

/**
 * @brief Print indexed complex single precision
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void ciprint(const int fold, const float_complex_indexed *X) {
  cmprint(fold, X, 2, X + fold * 2, 2);
}
