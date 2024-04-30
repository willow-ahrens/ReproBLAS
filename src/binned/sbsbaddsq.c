#include <binned.h>

/**
 * @brief Add binned single precision scaled sums of squares (Y += X)
 *
 * Performs the operation Y += X, where X and Y represent scaled sums of squares.
 *
 * @param fold the fold of the binned types
 * @param scaleX scale of X (scaleX == #binned_sscale (Z) for some @c float Z)
 * @param X binned scalar X
 * @param scaleY scale of Y (scaleY == #binned_sscale (Z) for some @c float Z)
 * @param Y binned scalar Y
 *
 * @return updated scale of Y
 *
 * @author Willow Ahrens
 * @date   2 Dec 2015
 */
float binned_sbsbaddsq(const int fold, const float scaleX, const float_binned *X, const float scaleY, float_binned *Y) {
return binned_smsmaddsq(fold, scaleX, X, 1, X + fold, 1, scaleY, Y, 1, Y + fold, 1);
}
