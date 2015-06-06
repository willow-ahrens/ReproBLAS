#include <indexed.h>

/**
 * @brief indexed double precision allocation
 *
 * @param fold the fold of the indexed type
 * @return a freshly allocated indexed type. (free with @c free())
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
double_indexed *dialloc(const int fold){
  return (double_indexed*)malloc(disize(fold));
}
