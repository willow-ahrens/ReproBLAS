#include <binned.h>

/**
 * @brief Convert binned double precision to double precision (X -> Y)
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 * @return scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
double binned_ddbconv(const int fold, const double_binned *X) {
  return binned_ddmconv(fold, X, 1, X + fold, 1);
}
