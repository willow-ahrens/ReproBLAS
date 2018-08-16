#include <binned.h>

/**
 * @brief binned single precision size
 *
 * @param fold the fold of the binned type
 * @return the size (in bytes) of the binned type
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
size_t binned_sbsbze(const int fold){
  return 2*fold*sizeof(float);
}
