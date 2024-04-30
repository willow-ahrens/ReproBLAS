#include <binned.h>

/**
 * @brief Set binned complex double precision to binned double precision (Y = X)
 *
 * Performs the operation Y = X
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_zbdbset(const int fold, const double_binned *X, double_complex_binned *Y){
  binned_zmdmset(fold, X, 1, X + fold, 1, Y, 1, Y + 2 * fold, 1);
}
