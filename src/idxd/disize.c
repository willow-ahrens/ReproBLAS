#include <idxd.h>

/**
 * @brief indexed double precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in bytes) of the indexed type
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
size_t idxd_disize(const int fold){
  return 2*fold*sizeof(double);
}
