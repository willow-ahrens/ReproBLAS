#include <binned.h>

/**
 * @internal
 * @brief Update manually specified binned complex single precision with complex single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value of real and imaginary components less than absolute value of real and imaginary components of X respectively.
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
void binned_cmcupdate(const int fold, const void *X, float* priY, const int incpriY, float* carY, const int inccarY) {
  binned_smsupdate(fold, ((float*)X)[0], priY, 2 * incpriY, carY, 2 * inccarY);
  binned_smsupdate(fold, ((float*)X)[1], priY + 1, 2 * incpriY, carY + 1, 2 * inccarY);
}
