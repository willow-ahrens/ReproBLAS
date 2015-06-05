#include <indexed.h>

/**
 * @internal
 * @brief Set manually specified indexed complex single precision (Y = X)
 *
 * Performs the operation Y = X
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmcmset(const int fold, const float *manX, const int incmanX, const float *carX, const int inccarX, float *manY, const int incmanY, float *carY, const int inccarY){
  smsmset(fold, manX, 2 * incmanX, carX, 2 * inccarX, manY, 2 * incmanY, carY, 2 * inccarY);
  smsmset(fold, manX + 1, 2 * incmanX, carX + 1, 2 * inccarX, manY + 1, 2 * incmanY, carY + 1, 2 * inccarY);
}
