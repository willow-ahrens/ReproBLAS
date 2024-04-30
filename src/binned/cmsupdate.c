#include <binned.h>

/**
 * @internal
 * @brief Update manually specified binned complex single precision with single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value less than X
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
void binned_cmsupdate(const int fold, const float X, float* priY, const int incpriY, float* carY, const int inccarY) {
  binned_smsupdate(fold, X, priY, 2 * incpriY, carY, 2 * inccarY);
  binned_smsupdate(fold, X, priY + 1, 2 * incpriY, carY + 1, 2 * inccarY);
}
