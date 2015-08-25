#include <string.h>

#include <idxd.h>

/**
 * @brief Set indexed complex single precision (Y = X)
 *
 * Performs the operation Y = X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void idxd_ciciset(const int fold, const float_complex_indexed *X, float_complex_indexed *Y){
  memcpy(Y, X, idxd_cisize(fold));
}
