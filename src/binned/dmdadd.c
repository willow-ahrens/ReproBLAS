#include <binned.h>

/**
 * @internal
 * @brief  Add double precision to manually specified binned double precision (Y += X)
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
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_dmdadd(const int fold, const double X, double *priY, const int incpriY, double *carY, const int inccarY){
  binned_dmdupdate(fold, X, priY, incpriY, carY, inccarY);
  binned_dmddeposit(fold, X, priY, incpriY);
  binned_dmrenorm(fold, priY, incpriY, carY, inccarY);
}
