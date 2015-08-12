#include <idxd.h>

/**
 * @brief Convert double precision to indexed double precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void didconv(const int fold, const double X, double_indexed *Y) {
  dmdconv(fold, X, Y, 1, Y + fold, 1);
}
