#include <binned.h>

/**
 * @brief Set binned complex single precision to binned single precision (Y = X)
 *
 * Performs the operation Y = X
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_cbsbset(const int fold, const float_binned *X, float_complex_binned *Y){
  binned_cmsmset(fold, X, 1, X + fold, 1, Y, 1, Y + 2 * fold, 1);
}
