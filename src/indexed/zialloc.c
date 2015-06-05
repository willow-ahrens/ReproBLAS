#include <stdlib.h>

#include <indexed.h>

/**
 * @brief indexed complex double precision allocation
 *
 * @param fold the fold of the indexed type
 * @return a freshly allocated indexed type. (free with @c free())
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
double_complex_indexed *zialloc(const int fold){
  return (double_complex_indexed*)malloc(zisize(fold));
}
