#include <binned.h>

/**
 * @brief Renormalize binned complex double precision
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
void binned_zbrenorm(const int fold, double_complex_binned *X) {
  binned_zmrenorm(fold, X, 1, X + 2 * fold, 1);
}
