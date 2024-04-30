#include <binned.h>

/**
 * @internal
 * @brief  Add complex single precision to manually specified binned complex single precision (Y += X)
 *
 * Performs the operation Y += X on an binned type Y
 *
 * @param fold the fold of the binned types
 * @param X scalar X
 * @param priY Y's primary vector
 * @param incpriY stride within Y's primary vector (use every incpriY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_cmcadd(const int fold, const void *X, float *priY, const int incpriY, float *carY, const int inccarY){
  binned_cmcupdate(fold, X, priY, incpriY, carY, inccarY);
  binned_cmcdeposit(fold, X, priY, incpriY);
  binned_cmrenorm(fold, priY, incpriY, carY, inccarY);
}
