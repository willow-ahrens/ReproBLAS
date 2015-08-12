#include <idxd.h>

/**
 * @brief  Add single precision to indexed single precision (Y += X)
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
void idxd_sisadd(const int fold, const float X, float_indexed *Y){
  idxd_smsadd(fold, X, Y, 1, Y + fold, 1);
}
