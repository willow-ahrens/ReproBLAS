#include <binned.h>

/**
 * @internal
 * @brief Set manually specified binned complex double precision to manually specified binned double precision (Y = X)
 *
 * Performs the operation Y = X
 *
 * @param fold the fold of the binned types
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param priY Y's primary vector
 * @param incpriY stride within Y's primary vector (use every incpriY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   27 Apr 2015
 */
void binned_zmdmset(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, double *priY, const int incpriY, double *carY, const int inccarY){
  binned_dmdmset(fold, priX, incpriX, carX, inccarX, priY, 2 * incpriY, carY, 2 * inccarY);
  binned_dmsetzero(fold, priY + 1, 2 * incpriY, carY + 1, 2 * inccarY);
}
