#include <binned.h>

/**
 * @brief  Add single precision to binned single precision (Y += X)
 *
 * Performs the operation Y += X on an binned type Y
 *
 * @param fold the fold of the binned types
 * @param X scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_sbsadd(const int fold, const float X, float_binned *Y){
  binned_smsadd(fold, X, Y, 1, Y + fold, 1);
}
