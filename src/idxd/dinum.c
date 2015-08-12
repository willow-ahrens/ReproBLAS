#include <idxd.h>

/**
 * @brief indexed double precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in @c double) of the indexed type
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int dinum(const int fold){
  return 2*fold;
}
