#include <binned.h>

/**
 * @brief binned complex double precision allocation
 *
 * @param fold the fold of the binned type
 * @return a freshly allocated binned type. (free with @c free())
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
double_complex_binned *binned_zballoc(const int fold){
  return (double_complex_binned*)malloc(binned_zbsize(fold));
}
