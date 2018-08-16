#include <binned.h>

/**
 * @brief  Add complex double precision to binned complex double precision (Y += X)
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
void binned_zbzadd(const int fold, const void *X, double_complex_binned *Y){
  binned_zmzadd(fold, X, Y, 1, Y + 2 * fold, 1);
}
