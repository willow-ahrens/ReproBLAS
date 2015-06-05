#include <indexed.h>

/**
 * @internal
 * @brief Convert manually specified indexed complex single precision to complex single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param conv scalar return
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void ccmconv_sub(const int fold, const float *manX, const int incmanX, const float *carX, const int inccarX, void *conv) {
  ((float*)conv)[0] = ssmconv(fold, manX, incmanX * 2, carX, inccarX + 1);
  ((float*)conv)[1] = ssmconv(fold, manX + 1, incmanX * 2, carX + 1, inccarX + 1);
}
