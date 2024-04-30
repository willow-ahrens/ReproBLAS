#include <binned.h>

/**
 * @brief  Add binned complex double precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_zbzbadd(const int fold, const double_complex_binned *X, double_complex_binned *Y){
  binned_zmzmadd(fold, X, 1, X + 2 * fold, 1, Y, 1, Y + 2 * fold, 1);
}
