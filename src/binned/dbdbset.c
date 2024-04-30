#include <string.h>

#include <binned.h>

/**
 * @brief Set binned double precision (Y = X)
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
void binned_dbdbset(const int fold, const double_binned *X, double_binned *Y){
  memcpy(Y, X, binned_dbsize(fold));
}
