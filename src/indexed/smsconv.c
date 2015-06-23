#include <indexed.h>

/**
 * @internal
 * @brief Convert single precision to manually specified indexed single precision (X -> Y)
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
 * @date   30 Apr 2015
 */
void smsconv(const int fold, const float X, float* priY, const int incpriY, float* carY, const int inccarY) {
  smsetzero(fold, priY, incpriY, carY, inccarY);
  smsadd(fold, X, priY, incpriY, carY, inccarY);
}
