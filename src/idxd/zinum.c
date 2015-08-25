#include <idxd.h>

/**
 * @brief indexed complex double precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in @c double) of the indexed type
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int idxd_zinum(const int fold){
  return 4*fold;
}
