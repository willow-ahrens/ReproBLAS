#include <binned.h>

/**
 * @brief Update binned complex single precision with single precision (X -> Y)
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
void binned_cbsupdate(const int fold, const float X, float_complex_binned *Y) {
  binned_cmsupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}
