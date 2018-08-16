#include <binned.h>

/**
 * @brief Negate binned double precision (X = -X)
 *
 * Performs the operation X = -X
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_dbnegate(const int fold, double_binned* X){
  binned_dmnegate(fold, X, 1, X + fold, 1);
}
