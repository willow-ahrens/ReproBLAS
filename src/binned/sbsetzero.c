#include <string.h>

#include <binned.h>

/**
 * @brief Set binned single precision to 0 (X = 0)
 *
 * Performs the operation X = 0
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_sbsetzero(const int fold, float_binned *X){
  memset(X, 0, binned_sbsbze(fold));
}
