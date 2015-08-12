#include <idxd.h>

/**
 * @brief indexed single precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in bytes) of the indexed type
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
size_t sisize(const int fold){
  return 2*fold*sizeof(float);
}
