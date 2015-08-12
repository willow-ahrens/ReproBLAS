#include <idxd.h>

/**
 * @brief  Add indexed double precision (Y += X)
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
void didiadd(const int fold, const double_indexed *X, double_indexed *Y){
  dmdmadd(fold, X, 1, X + fold, 1, Y, 1, Y + fold, 1);
}
