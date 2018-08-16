#include <binned.h>

/**
 * @brief Print binned double precision
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_dbprint(const int fold, const double_binned *X) {
  binned_dmprint(fold, X, 1, X + fold, 1);
}
