#include <binned.h>

/**
 * @brief  Add complex single precision to binned complex single precision (Y += X)
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
void binned_cbcadd(const int fold, const void *X, float_complex_binned *Y){
  binned_cmcadd(fold, X, Y, 1, Y + 2 * fold, 1);
}
