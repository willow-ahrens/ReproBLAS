#include <binned.h>

/**
 * @brief Convert single precision to binned single precision (X -> Y)
 *
 * @param fold the fold of the binned types
 * @param X scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_sbsconv(const int fold, const float X, float_binned *Y) {
  binned_smsconv(fold, X, Y, 1, Y + fold, 1);
}
