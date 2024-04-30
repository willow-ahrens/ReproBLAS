#include <binned.h>

/**
 * @brief binned single precision size
 *
 * @param fold the fold of the binned type
 * @return the size (in @c float) of the binned type
 *
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
int binned_sbnum(const int fold){
  return 2*fold;
}
