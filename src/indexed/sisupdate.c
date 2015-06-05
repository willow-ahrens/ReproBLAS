#include <indexed.h>

/**
 * @brief Update indexed single precision with single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value less than X
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void sisupdate(const int fold, const float X, float_indexed *Y) {
  smsupdate(fold, X, Y, 1, Y + fold, 1);
}
