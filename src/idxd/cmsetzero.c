#include <idxd.h>

/**
 * @internal
 * @brief Set manually specified indexed complex single precision to 0 (X = 0)
 *
 * Performs the operation X = 0
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
void cmsetzero(const int fold, float *priX, const int incpriX, float *carX, const int inccarX){
  smsetzero(fold, priX, 2 * incpriX, carX, 2 * inccarX);
  smsetzero(fold, priX + 1, 2 * incpriX, carX + 1, 2 * inccarX);
}
