#include <indexed.h>

/**
 * @brief Convert indexed complex double precision to complex double precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param conv scalar return
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zziconv_sub(const int fold, const double_complex_indexed *X, void *conv) {
  zzmconv_sub(fold, X, 1, X + 2 * fold, 1, conv);
}
