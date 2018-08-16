#include <binned.h>

/**
 * @brief Convert binned complex double precision to complex double precision (X -> Y)
 *
 * @param fold the fold of the binned types
 * @param X binned scalar X
 * @param conv scalar return
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_zzbconv_sub(const int fold, const double_complex_binned *X, void *conv) {
  binned_zzmconv_sub(fold, X, 1, X + 2 * fold, 1, conv);
}
