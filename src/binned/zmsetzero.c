#include <binned.h>

/**
 * @internal
 * @brief Set manually specified binned complex double precision to 0 (X = 0)
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
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void binned_zmsetzero(const int fold, double *priX, const int incpriX, double *carX, const int inccarX){
  binned_dmsetzero(fold, priX, 2 * incpriX, carX, 2 * inccarX);
  binned_dmsetzero(fold, priX + 1, 2 * incpriX, carX + 1, 2 * inccarX);
}
