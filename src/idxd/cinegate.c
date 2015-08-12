#include <idxd.h>

/**
 * @brief Negate indexed complex single precision (X = -X)
 *
 * Performs the operation X = -X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cinegate(const int fold, float_complex_indexed* X){
  cmnegate(fold, X, 1, X + 2 * fold, 1);
}
