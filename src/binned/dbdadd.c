#include <binned.h>

/**
 * @brief  Add double precision to binned double precision (Y += X)
 *
 * Performs the operation Y += X on an binned type Y
 *
 * @param fold the fold of the binned types
 * @param X scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_dbdadd(const int fold, const double X, double_binned *Y){
  binned_dmdadd(fold, X, Y, 1, Y + fold, 1);
}
