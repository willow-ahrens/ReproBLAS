#include <string.h>

#include <idxd.h>

/**
 * @brief Set indexed complex double precision (Y = X)
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
void ziziset(const int fold, const double_complex_indexed *X, double_complex_indexed *Y){
  memcpy(Y, X, zisize(fold));
}
