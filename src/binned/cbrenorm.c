#include <binned.h>

/**
 * @brief Renormalize binned complex single precision
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
void binned_cbrenorm(const int fold, float_complex_binned *X) {
  binned_cmrenorm(fold, X, 1, X + 2 * fold, 1);
}
