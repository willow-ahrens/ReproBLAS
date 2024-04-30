#include <binned.h>

/**
 * @brief Update binned single precision with single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value less than X
 *
 * @param fold the fold of the binned types
 * @param X scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_sbsupdate(const int fold, const float X, float_binned *Y) {
  binned_smsupdate(fold, X, Y, 1, Y + fold, 1);
}
