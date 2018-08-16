#include <binned.h>

/**
 * @brief Print binned single precision
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_sbprint(const int fold, const float_binned *X) {
  binned_smprint(fold, X, 1, X + fold, 1);
}
