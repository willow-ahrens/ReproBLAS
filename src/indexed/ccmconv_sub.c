#include <indexed.h>

/**
 * @internal
 * @brief Convert manually specified indexed complex single precision to complex single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param conv scalar return
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void ccmconv_sub(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, void *conv) {
  ((float*)conv)[0] = ssmconv(fold, priX, incpriX * 2, carX, inccarX + 1);
  ((float*)conv)[1] = ssmconv(fold, priX + 1, incpriX * 2, carX + 1, inccarX + 1);
}
