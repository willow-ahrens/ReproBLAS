#include <indexed.h>

/**
 * @internal
 * @brief  Add complex single precision to manually specified indexed complex single precision (Y += X)
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
void cmcadd(const int fold, const void *X, float *manY, const int incmanY, float *carY, const int inccarY){
  cmcupdate(fold, X, manY, incmanY, carY, inccarY);
  cmcdeposit(fold, X, manY, incmanY);
  cmrenorm(fold, manY, incmanY, carY, inccarY);
}
