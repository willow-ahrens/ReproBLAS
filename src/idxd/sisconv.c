#include <idxd.h>

/**
 * @brief Convert single precision to indexed single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void idxd_sisconv(const int fold, const float X, float_indexed *Y) {
  idxd_smsconv(fold, X, Y, 1, Y + fold, 1);
}
