#include <binned.h>

/**
 * @brief  Add binned single precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_sbsbadd(const int fold, const float_binned *X, float_binned *Y){
  binned_smsmadd(fold, X, 1, X + fold, 1, Y, 1, Y + fold, 1);
}
