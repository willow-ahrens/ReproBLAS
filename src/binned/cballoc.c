#include <binned.h>

/**
 * @brief binned complex single precision allocation
 *
 * @param fold the fold of the binned type
 * @return a freshly allocated binned type. (free with @c free())
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
float_complex_binned *binned_cballoc(const int fold){
  return (float_complex_binned*)malloc(binned_cbsize(fold));
}
