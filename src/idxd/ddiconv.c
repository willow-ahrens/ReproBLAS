#include <idxd.h>

/**
 * @brief Convert indexed double precision to double precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @return scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
double idxd_ddiconv(const int fold, const double_indexed *X) {
  return idxd_ddmconv(fold, X, 1, X + fold, 1);
}
