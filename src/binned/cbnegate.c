#include <binned.h>

/**
 * @brief Negate binned complex single precision (X = -X)
 *
 * Performs the operation X = -X
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_cbnegate(const int fold, float_complex_binned* X){
  binned_cmnegate(fold, X, 1, X + 2 * fold, 1);
}
