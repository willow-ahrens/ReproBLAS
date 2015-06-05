#include <indexed.h>

/**
 * @brief  Add complex single precision to indexed complex single precision (Y += X)
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
void cicadd(const int fold, const void *X, float_complex_indexed *Y){
  cmcadd(fold, X, Y, 1, Y + 2 * fold, 1);
}
