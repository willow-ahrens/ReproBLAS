#include <idxd.h>

/**
 * @brief  Add double precision to indexed double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void idxd_didadd(const int fold, const double X, double_indexed *Y){
  idxd_dmdadd(fold, X, Y, 1, Y + fold, 1);
}
