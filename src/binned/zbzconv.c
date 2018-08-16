#include <binned.h>

/**
 * @brief Convert complex double precision to binned complex double precision (X -> Y)
 *
 * @param fold the fold of the binned types
 * @param X scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_zbzconv(const int fold, const void *X, double_complex_binned *Y) {
  binned_zmzconv(fold, X, Y, 1, Y + 2 * fold, 1);
}
