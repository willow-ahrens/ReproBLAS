#include <idxd.h>

/**
 * @brief indexed single precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in @c float) of the indexed type
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int idxd_sinum(const int fold){
  return 2*fold;
}
