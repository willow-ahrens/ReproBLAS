#include <idxd.h>

/**
 * @brief Negate indexed double precision (X = -X)
 *
 * Performs the operation X = -X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void idxd_dinegate(const int fold, double_indexed* X){
  idxd_dmnegate(fold, X, 1, X + fold, 1);
}
