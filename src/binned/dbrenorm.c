#include <binned.h>

/**
 * @brief Renormalize binned double precision
 *
 * Renormalization keeps the primary vector within the necessary bins by shifting over to the carry vector
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_dbrenorm(const int fold, double_binned *X) {
  binned_dmrenorm(fold, X, 1, X + fold, 1);
}
