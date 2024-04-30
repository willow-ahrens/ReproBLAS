#include <binned.h>

/**
 * @internal
 * @brief Convert manually specified binned complex double precision to complex double precision (X -> Y)
 *
 * @param fold the fold of the binned types
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param conv scalar return
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_zzmconv_sub(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, void *conv) {
  ((double*)conv)[0] = binned_ddmconv(fold, priX, incpriX * 2, carX, inccarX + 1);
  ((double*)conv)[1] = binned_ddmconv(fold, priX + 1, incpriX * 2, carX + 1, inccarX + 1);
}
