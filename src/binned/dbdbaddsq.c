#include <binned.h>

/**
 * @brief Add binned double precision scaled sums of squares (Y += X)
 *
 * Performs the operation Y += X, where X and Y represent scaled sums of squares.
 *
 * @param fold the fold of the binned types
 * @param scaleX scale of X (scaleX == #binned_dscale (Z) for some @c double Z)
 * @param X binned scalar X
 * @param scaleY scale of Y (scaleY == #binned_dscale (Z) for some @c double Z)
 * @param Y binned scalar Y
 *
 * @return updated scale of Y
 *
 * @author Peter Ahrens
 * @date   2 Dec 2015
 */
double binned_dbdbaddsq(const int fold, const double scaleX, const double_binned *X, const double scaleY, double_binned *Y) {
return binned_dmdmaddsq(fold, scaleX, X, 1, X + fold, 1, scaleY, Y, 1, Y + fold, 1);
}
