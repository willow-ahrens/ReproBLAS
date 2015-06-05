#include <indexed.h>

/**
 * @brief Update indexed complex double precision with double precision (X -> Y)
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
void zidupdate(const int fold, const double X, double_complex_indexed *Y) {
  zmdupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}
