#include <indexed.h>

/**
 * @brief Negate indexed complex double precision (X = -X)
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
void zinegate(const int fold, double_complex_indexed* X){
  zmnegate(fold, X, 1, X + 2 * fold, 1);
}
