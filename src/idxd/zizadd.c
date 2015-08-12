#include <idxd.h>

/**
 * @brief  Add complex double precision to indexed complex double precision (Y += X)
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
void idxd_zizadd(const int fold, const void *X, double_complex_indexed *Y){
  idxd_zmzadd(fold, X, Y, 1, Y + 2 * fold, 1);
}
