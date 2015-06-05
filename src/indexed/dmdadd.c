#include <indexed.h>

/**
 * @internal
 * @brief  Add double precision to manually specified indexed double precision (Y += X)
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
void dmdadd(const int fold, const double X, double *manY, const int incmanY, double *carY, const int inccarY){
  dmdupdate(fold, X, manY, incmanY, carY, inccarY);
  dmddeposit(fold, X, manY, incmanY);
  dmrenorm(fold, manY, incmanY, carY, inccarY);
}
