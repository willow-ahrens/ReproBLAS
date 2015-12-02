#include <idxd.h>

/**
 * @brief Add indexed double precision scaled sums of squares (Y += X)
 *
 * Performs the operation Y += X, where X and Y represent scaled sums of squares.
 *
 * @param fold the fold of the indexed types
 * @param scaleX scale of X (scaleX == #idxd_dscale(Z) for some @c double Z)
 * @param X indexed scalar X
 * @param scaleY scale of Y (scaleY == #idxd_dscale(Z) for some @c double Z)
 * @param Y indexed scalar Y
 *
 * @return updated scale of Y
 *
 * @author Peter Ahrens
 * @date   2 Dec 2015
 */
double idxd_didiaddsq(const int fold, const double scaleX, const double_indexed *X, const double scaleY, double_indexed *Y) {
return idxd_dmdmaddsq(fold, scaleX, X, 1, X + fold, 1, scaleY, Y, 1, Y + fold, 1);
}
