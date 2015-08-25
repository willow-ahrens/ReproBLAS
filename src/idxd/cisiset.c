#include <idxd.h>

/**
 * @brief Set indexed complex single precision to indexed single precision (Y = X)
 *
 * Performs the operation Y = X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void idxd_cisiset(const int fold, const float_indexed *X, float_complex_indexed *Y){
  idxd_cmsmset(fold, X, 1, X + fold, 1, Y, 1, Y + 2 * fold, 1);
}
