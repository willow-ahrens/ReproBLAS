#include <string.h>

#include <binned.h>

/**
 * @brief Set binned complex single precision (Y = X)
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
void binned_cbcbset(const int fold, const float_complex_binned *X, float_complex_binned *Y){
  memcpy(Y, X, binned_cbsize(fold));
}
