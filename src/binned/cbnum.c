#include <binned.h>

/**
 * @brief binned complex single precision size
 *
 * @param fold the fold of the binned type
 * @return the size (in @c float) of the binned type
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int binned_cbnum(const int fold){
  return 4*fold;
}
