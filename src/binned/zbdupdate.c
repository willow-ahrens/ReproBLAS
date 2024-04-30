#include <binned.h>

/**
 * @brief Update binned complex double precision with double precision (X -> Y)
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
void binned_zbdupdate(const int fold, const double X, double_complex_binned *Y) {
  binned_zmdupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}
