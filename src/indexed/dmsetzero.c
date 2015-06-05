#include <indexed.h>

/**
 * @internal
 * @brief Set manually specified indexed double precision to 0 (X = 0)
 *
 * Performs the operation X = 0
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dmsetzero(const int fold, double *manX, const int incmanX, double *carX, const int inccarX){
  int i;
  for(i = 0; i < fold; i++){
    manX[i * incmanX] = 0.0;
    carX[i * inccarX] = 0.0;
  }
}
