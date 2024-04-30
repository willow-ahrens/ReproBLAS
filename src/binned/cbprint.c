#include <binned.h>

/**
 * @brief Print binned complex single precision
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_cbprint(const int fold, const float_complex_binned *X) {
  binned_cmprint(fold, X, 2, X + fold * 2, 2);
}
