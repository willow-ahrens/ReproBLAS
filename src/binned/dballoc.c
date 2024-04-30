#include <binned.h>

/**
 * @brief binned double precision allocation
 *
 * @param fold the fold of the binned type
 * @return a freshly allocated binned type. (free with @c free())
 *
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
double_binned *binned_dballoc(const int fold){
  return (double_binned*)malloc(binned_dbsize(fold));
}
