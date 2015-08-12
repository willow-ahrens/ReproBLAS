#include <idxd.h>

/**
 * @brief Convert complex single precision to indexed complex single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cicconv(const int fold, const void *X, float_complex_indexed *Y) {
  cmcconv(fold, X, Y, 1, Y + 2 * fold, 1);
}
