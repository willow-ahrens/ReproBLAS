#include <string.h>

#include <binned.h>

/**
 * @brief Set binned complex double precision (Y = X)
 *
 * Performs the operation Y = X
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_zbzbset(const int fold, const double_complex_binned *X, double_complex_binned *Y){
  memcpy(Y, X, binned_zbsize(fold));
}
