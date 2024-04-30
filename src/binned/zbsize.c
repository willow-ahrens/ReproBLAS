#include <binned.h>

/**
 * @brief binned complex double precision size
 *
 * @param fold the fold of the binned type
 * @return the size (in bytes) of the binned type
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
size_t binned_zbsize(const int fold){
  return 4*fold*sizeof(double);
}
