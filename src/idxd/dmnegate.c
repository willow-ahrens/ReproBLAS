#include <idxd.h>

/**
 * @internal
 * @brief Negate manually specified indexed double precision (X = -X)
 *
 * Performs the operation X = -X
 *
 * @param fold the fold of the indexed types
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dmnegate(const int fold, double* priX, const int incpriX, double* carX, const int inccarX) {
  int i;
  const double *bins;
  if(priX[0] != 0.0){
    bins = dmbins(dmindex(priX));
    for (i = 0; i < fold; i++) {
      priX[i * incpriX] = bins[i] - (priX[i * incpriX] - bins[i]);
      carX[i * inccarX] = -carX[i * inccarX];
    }
  }
}
