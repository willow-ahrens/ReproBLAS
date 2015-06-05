#include <indexed.h>

/**
 * @internal
 * @brief Update manually specified indexed complex single precision with single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value less than X
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
void cmsupdate(const int fold, const float X, float* manY, const int incmanY, float* carY, const int inccarY) {
  smsupdate(fold, X, manY, 2 * incmanY, carY, 2 * inccarY);
  smsupdate(fold, X, manY + 1, 2 * incmanY, carY + 1, 2 * inccarY);
}
