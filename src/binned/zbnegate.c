#include <binned.h>

/**
 * @brief Negate binned complex double precision (X = -X)
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
void binned_zbnegate(const int fold, double_complex_binned* X){
  binned_zmnegate(fold, X, 1, X + 2 * fold, 1);
}
