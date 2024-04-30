#include <string.h>

#include <binned.h>

/**
 * @brief Set binned single precision (Y = X)
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
void binned_sbsbset(const int fold, const float_binned *X, float_binned *Y){
  memcpy(Y, X, binned_sbsbze(fold));
}
