#include <idxd.h>

/**
 * @brief  Add indexed single precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void sisiadd(const int fold, const float_indexed *X, float_indexed *Y){
  smsmadd(fold, X, 1, X + fold, 1, Y, 1, Y + fold, 1);
}
