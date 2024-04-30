#include <binned.h>

/**
 * @brief Update binned complex double precision with complex double precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value of real and imaginary components less than absolute value of real and imaginary components of X respectively.
 *
 * @param fold the fold of the binned types
 * @param X scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_zbzupdate(const int fold, const void *X, double_complex_binned *Y) {
  binned_zmzupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}
