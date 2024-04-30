#include <binned.h>

/**
 * @brief Renormalize binned single precision
 *
 * Renormalization keeps the primary vector within the necessary bins by shifting over to the carry vector
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_sbrenorm(const int fold, float_binned *X) {
  binned_smrenorm(fold, X, 1, X + fold, 1);
}
