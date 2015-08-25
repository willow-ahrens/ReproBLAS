#include <idxd.h>

/**
 * @brief Update indexed complex single precision with single precision (X -> Y)
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
void idxd_cisupdate(const int fold, const float X, float_complex_indexed *Y) {
  idxd_cmsupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}
