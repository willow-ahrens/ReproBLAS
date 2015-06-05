#include <indexed.h>

/**
 * @internal
 * @brief  Add single precision to manually specified indexed single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
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
void smsadd(const int fold, const float X, float *manY, const int incmanY, float *carY, const int inccarY){
  smsupdate(fold, X, manY, incmanY, carY, inccarY);
  smsdeposit(fold, X, manY, incmanY);
  smrenorm(fold, manY, incmanY, carY, inccarY);
}
