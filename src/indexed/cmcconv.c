#include <indexed.h>

/**
 * @internal
 * @brief Convert complex single precision to manually specified indexed complex single precision (X -> Y)
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
void cmcconv(const int fold, const void *X, float *manY, const int incmanY, float *carY, const int inccarY) {
  smsconv(fold, ((float*)X)[0], manY, incmanY * 2, carY, inccarY * 2);
  smsconv(fold, ((float*)X)[1], manY + 1, incmanY * 2, carY + 1, inccarY * 2);
}
