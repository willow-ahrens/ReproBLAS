#include <binned.h>

/**
 * @internal
 * @brief Convert complex double precision to manually specified binned complex double precision (X -> Y)
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
void binned_zmzconv(const int fold, const void *X, double *priY, const int incpriY, double *carY, const int inccarY) {
  binned_dmdconv(fold, ((double*)X)[0], priY, incpriY * 2, carY, inccarY * 2);
  binned_dmdconv(fold, ((double*)X)[1], priY + 1, incpriY * 2, carY + 1, inccarY * 2);
}
