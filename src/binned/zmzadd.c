#include <binned.h>

/**
 * @internal
 * @brief  Add complex double precision to manually specified binned complex double precision (Y += X)
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
void binned_zmzadd(const int fold, const void *X, double *priY, const int incpriY, double *carY, const int inccarY){
  binned_zmzupdate(fold, X, priY, incpriY, carY, inccarY);
  binned_zmzdeposit(fold, X, priY, incpriY);
  binned_zmrenorm(fold, priY, incpriY, carY, inccarY);
}
