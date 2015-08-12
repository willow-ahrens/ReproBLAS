#include <idxd.h>

/**
 * @brief indexed complex single precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in @c float) of the indexed type
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int idxd_cinum(const int fold){
  return 4*fold;
}
