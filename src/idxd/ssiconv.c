#include <idxd.h>

/**
 * @brief Convert indexed single precision to single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @return scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
float ssiconv(const int fold, const float_indexed *X) {
  return ssmconv(fold, X, 1, X + fold, 1);
}
