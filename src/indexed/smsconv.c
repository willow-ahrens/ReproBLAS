#include <indexed.h>

/**
 * @internal
 * @brief Convert single precision to manually specified indexed single precision (X -> Y)
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
 * @date   30 Apr 2015
 */
void smsconv(const int fold, const float X, float* manY, const int incmanY, float* carY, const int inccarY) {
  smsetzero(fold, manY, incmanY, carY, inccarY);
  smsadd(fold, X, manY, incmanY, carY, inccarY);
}
