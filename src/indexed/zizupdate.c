#include <indexed.h>

/**
 * @brief Update indexed complex double precision with complex double precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value of real and imaginary components less than absolute value of real and imaginary components of X respectively.
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zizupdate(const int fold, const void *X, double_complex_indexed *Y) {
  zmzupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}
