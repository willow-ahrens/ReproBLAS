#include <indexed.h>

/**
 * @internal
 * @brief Convert complex double precision to manually specified indexed complex double precision (X -> Y)
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
void zmzconv(const int fold, const void *X, double *manY, const int incmanY, double *carY, const int inccarY) {
  dmdconv(fold, ((double*)X)[0], manY, incmanY * 2, carY, inccarY * 2);
  dmdconv(fold, ((double*)X)[1], manY + 1, incmanY * 2, carY + 1, inccarY * 2);
}
