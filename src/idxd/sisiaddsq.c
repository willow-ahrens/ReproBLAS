#include <idxd.h>

/**
 * @brief Add indexed single precision scaled sums of squares (Y += X)
 *
 * Performs the operation Y += X, where X and Y represent scaled sums of squares.
 *
 * @param fold the fold of the indexed types
 * @param scaleX scale of X (scaleX == #idxd_sscale(Z) for some @c float Z)
 * @param X indexed scalar X
 * @param scaleY scale of Y (scaleY == #idxd_sscale(Z) for some @c float Z)
 * @param Y indexed scalar Y
 *
 * @return updated scale of Y
 *
 * @author Peter Ahrens
 * @date   2 Dec 2015
 */
float idxd_sisiaddsq(const int fold, const float scaleX, const float_indexed *X, const float scaleY, float_indexed *Y) {
return idxd_smsmaddsq(fold, scaleX, X, 1, X + fold, 1, scaleY, Y, 1, Y + fold, 1);
}
