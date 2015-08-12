#include <string.h>

#include <idxd.h>

/**
 * @brief Set indexed double precision (Y = X)
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
void didiset(const int fold, const double_indexed *X, double_indexed *Y){
  memcpy(Y, X, disize(fold));
}
