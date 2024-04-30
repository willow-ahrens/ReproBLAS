#include <math.h>

#include <binned.h>

#include "../common/common.h"

/**
 * @internal
 * @brief  Add manually specified binned double precision (Y += X)
 *
 * Performs the operation Y += X
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
void binned_dmdmadd(const int fold, const double *priX, const int incpriX, const double *carX, const int inccarX, double* priY, const int incpriY, double* carY, const int inccarY) {
  int i;
  int shift;
  int X_index;
  int Y_index;
  const double *bins;

  if (priX[0] == 0.0)
    return;

  if (priY[0] == 0.0) {
    for (i = 0; i < fold; i++) {
      priY[i*incpriY] = priX[i*incpriX];
      carY[i*inccarY] = carX[i*inccarX];
    }
    return;
  }

  if (ISNANINF(priX[0]) || ISNANINF(priY[0])){
    priY[0] += priX[0];
    return;
  }

  X_index = binned_dmindex(priX);
  Y_index = binned_dmindex(priY);
  shift = Y_index - X_index;
  if(shift > 0){
    bins = binned_dmbins(Y_index);
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
    bins = binned_dmbins(X_index);
    //shift X upwards and add X to Y
    for (i = 0 - shift; i < fold; i++) {
      priY[i*incpriY] += priX[(i + shift)*incpriX] - bins[i + shift];
      carY[i*inccarY] += carX[(i + shift)*inccarX];
    }
  }

  binned_dmrenorm(fold, priY, incpriY, carY, inccarY);
}
