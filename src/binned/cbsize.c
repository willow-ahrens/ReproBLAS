#include <binned.h>

/**
 * @brief binned complex single precision size
 *
 * @param fold the fold of the binned type
 * @return the size (in bytes) of the binned type
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
size_t binned_cbsize(const int fold){
  return 4*fold*sizeof(float);
}
