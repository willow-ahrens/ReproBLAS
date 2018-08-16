#include <binned.h>

/**
 * @brief Convert binned single precision to single precision (X -> Y)
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 * @return scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
float binned_ssbconv(const int fold, const float_binned *X) {
  return binned_ssmconv(fold, X, 1, X + fold, 1);
}
