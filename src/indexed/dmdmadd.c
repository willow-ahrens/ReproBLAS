#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief  Add manually specified indexed double precision (Y += X)
 *
 * Performs the operation Y += X
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
void dmdmadd(const int fold, const double *manX, const int incmanX, const double *carX, const int inccarX, double* manY, const int incmanY, double* carY, const int inccarY) {
  int i;
  int shift;
  int X_index;
  int Y_index;
  const double *bins;

  if (manX[0] == 0.0)
    return;

  if (manY[0] == 0.0) {
    for (i = 0; i < fold; i++) {
      manY[i*incmanY] = manX[i*incmanX];
      carY[i*inccarY] = carX[i*inccarX];
    }
    return;
  }

  if (ISNANINF(manX[0]) || ISNANINF(manY[0])){
    manY[0] += manX[0];
    return;
  }

  X_index = dmindex(manX);
  Y_index = dmindex(manY);
  shift = Y_index - X_index;
  if(shift > 0){
    bins = dmbins(Y_index);
    //shift Y upwards and add X to Y
    for (i = fold - 1; i >= shift; i--) {
      manY[i*incmanY] = manX[i*incmanX] + (manY[(i - shift)*incmanY] - bins[i - shift]);
      carY[i*inccarY] = carX[i*inccarX] + carY[(i - shift)*inccarY];
    }
    for (i = 0; i < shift && i < fold; i++) {
      manY[i*incmanY] = manX[i*incmanX];
      carY[i*inccarY] = carX[i*inccarX];
    }
  }else{
    bins = dmbins(X_index);
    //shift X upwards and add X to Y
    for (i = 0 - shift; i < fold; i++) {
      manY[i*incmanY] += manX[(i + shift)*incmanX] - bins[i + shift];
      carY[i*inccarY] += carX[(i + shift)*inccarX];
    }
  }

  dmrenorm(fold, manY, incmanY, carY, inccarY);
}
