#include <binned.h>

/**
 * @brief binned single precision allocation
 *
 * @param fold the fold of the binned type
 * @return a freshly allocated binned type. (free with @c free())
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
float_binned *binned_sballoc(const int fold){
  return (float_binned*)malloc(binned_sbsbze(fold));
}
