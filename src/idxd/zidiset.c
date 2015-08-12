#include <idxd.h>

/**
 * @brief Set indexed complex double precision to indexed double precision (Y = X)
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
void zidiset(const int fold, const double_indexed *X, double_complex_indexed *Y){
  zmdmset(fold, X, 1, X + fold, 1, Y, 1, Y + 2 * fold, 1);
}
