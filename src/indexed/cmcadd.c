#include <indexed.h>

/**
 * @internal
 * @brief  Add complex single precision to manually specified indexed complex single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param priY Y's primary vector
 * @param incpriY stride within Y's primary vector (use every incpriY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmcadd(const int fold, const void *X, float *priY, const int incpriY, float *carY, const int inccarY){
  cmcupdate(fold, X, priY, incpriY, carY, inccarY);
  cmcdeposit(fold, X, priY, incpriY);
  cmrenorm(fold, priY, incpriY, carY, inccarY);
}
