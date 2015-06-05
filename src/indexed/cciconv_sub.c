#include <indexed.h>

/**
 * @brief Convert indexed complex single precision to complex single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param conv scalar return
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cciconv_sub(const int fold, const float_complex_indexed *X, void *conv) {
  ccmconv_sub(fold, X, 1, X + 2 * fold, 1, conv);
}
