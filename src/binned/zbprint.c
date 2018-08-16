#include <binned.h>

/**
 * @brief Print binned complex double precision
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_zbprint(const int fold, const double_complex_binned *X){
  binned_zmprint(fold, X, 2, X + 2 * fold, 2);
}
