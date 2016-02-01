#include <idxd.h>

/**
 * @brief Print indexed complex double precision
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void idxd_ziprint(const int fold, const double_complex_indexed *X){
  idxd_zmprint(fold, X, 2, X + 2 * fold, 2);
}
