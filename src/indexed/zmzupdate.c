#include <indexed.h>

/**
 * @internal
 * @brief Update manually specified indexed complex double precision with complex double precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value of real and imaginary components less than absolute value of real and imaginary components of X respectively.
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zmzupdate(const int fold, const void *X, double* manY, const int incmanY, double* carY, const int inccarY) {
  dmdupdate(fold, ((double*)X)[0], manY, 2 * incmanY, carY, 2 * inccarY);
  dmdupdate(fold, ((double*)X)[1], manY + 1, 2 * incmanY, carY + 1, 2 * inccarY);
}
