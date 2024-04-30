#include <binned.h>

/**
 * @brief Convert double precision to binned double precision (X -> Y)
 *
 * @param fold the fold of the binned types
 * @param X scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_dbdconv(const int fold, const double X, double_binned *Y) {
  binned_dmdconv(fold, X, Y, 1, Y + fold, 1);
}
