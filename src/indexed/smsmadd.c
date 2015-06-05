#include <math.h>

#include <indexed.h>

/**
 * @internal
 * @brief  Add manually specified indexed single precision (Y += X)
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
void smsmadd(const int fold, const float *manX, const int incmanX, const float *carX, const int inccarX, float* manY, const int incmanY, float* carY, const int inccarY) {
  int i;
  int shift;

  if (manX[0] == 0.0)
    return;

  if (manY[0] == 0.0) {
    for (i = 0; i < fold; i++) {
      manY[i*incmanY] = manX[i*incmanX];
      carY[i*inccarY] = carX[i*inccarX];
    }
    return;
  }

  if (isinf(manX[0]) || isnan(manX[0]) || isinf(manY[0]) || isnan(manY[0])) {
    manY[0] += manX[0];
    return;
  }

  shift = smindex(manY) - smindex(manX);
  if(shift > 0){
    //shift Y upwards and add X to Y
    for (i = fold - 1; i >= shift; i--) {
      manY[i*incmanY] = manX[i*incmanX] + (manY[(i - shift)*incmanY] - 1.5*ufpf(manY[(i - shift)*incmanY]));
      carY[i*inccarY] = carX[i*inccarX] + carY[(i - shift)*inccarY];
    }
    for (i = 0; i < shift && i < fold; i++) {
      manY[i*incmanY] = manX[i*incmanX];
      carY[i*inccarY] = carX[i*inccarX];
    }
  }else{
    //shift X upwards and add X to Y
    for (i = 0 - shift; i < fold; i++) {
      manY[i*incmanY] += manX[(i + shift)*incmanX] - 1.5*ufpf(manX[(i + shift)*incmanX]);
      carY[i*inccarY] += carX[(i + shift)*inccarX];
    }
  }

  smrenorm(fold, manY, incmanY, carY, inccarY);
}
