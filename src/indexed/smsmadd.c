#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief  Add manually specified indexed single precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
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
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void smsmadd(const int fold, const float *priX, const int incpriX, const float *carX, const int inccarX, float* priY, const int incpriY, float* carY, const int inccarY) {
  int i;
  int shift;
  int X_index;
  int Y_index;
  const float *bins;

  if (priX[0] == 0.0)
    return;

  if (priY[0] == 0.0) {
    for (i = 0; i < fold; i++) {
      priY[i*incpriY] = priX[i*incpriX];
      carY[i*inccarY] = carX[i*inccarX];
    }
    return;
  }

  if (ISNANINFF(priX[0]) || ISNANINFF(priY[0])){
    priY[0] += priX[0];
    return;
  }

  X_index = smindex(priX);
  Y_index = smindex(priY);
  shift = Y_index - X_index;
  if(shift > 0){
    bins = smbins(Y_index);
    //shift Y upwards and add X to Y
    for (i = fold - 1; i >= shift; i--) {
      priY[i*incpriY] = priX[i*incpriX] + (priY[(i - shift)*incpriY] - bins[i - shift]);
      carY[i*inccarY] = carX[i*inccarX] + carY[(i - shift)*inccarY];
    }
    for (i = 0; i < shift && i < fold; i++) {
      priY[i*incpriY] = priX[i*incpriX];
      carY[i*inccarY] = carX[i*inccarX];
    }
  }else{
    bins = smbins(X_index);
    //shift X upwards and add X to Y
    for (i = 0 - shift; i < fold; i++) {
      priY[i*incpriY] += priX[(i + shift)*incpriX] - bins[i + shift];
      carY[i*inccarY] += carX[(i + shift)*inccarX];
    }
  }

  smrenorm(fold, priY, incpriY, carY, inccarY);
}
