#include <binned.h>

/**
 * @internal
 * @brief Set manually specified binned double precision to 0 (X = 0)
 *
 * Performs the operation X = 0
 *
 * @param fold the fold of the binned types
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_dmsetzero(const int fold, double *priX, const int incpriX, double *carX, const int inccarX){
  int i;
  for(i = 0; i < fold; i++){
    priX[i * incpriX] = 0.0;
    carX[i * inccarX] = 0.0;
  }
}
