/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

void dmdupdate(double X, double* repY, int increpY, double* carY, int inccarY, int fold) {
  if (X == 0 || isnan(repY[0]) || isinf(repY[0]))
    return;

/*
  if (repY[0] == 0.0) {
    dmbound(dindex(X), repY, increpY, fold);
    for (int i = 0; i < fold; i++) {
      carY[i * inccarY] = 0.0;
    }
    return;
  }
*/

  int X_index = dindex(X);
  int d = diindex(repY) - X_index;
  if(d > 0){
    for(int i = fold - 1; i >= d; i--){
      repY[i * increpY] = repY[(i - d) * increpY];
      carY[i * inccarY] = carY[(i - d) * inccarY];
    }
    dmbound(X_index, repY, increpY, MIN(d, fold));
    for(int i = 0; i < d && i < fold; i++){
      carY[i * inccarY] = 0.0;
    }
  }
}

void didupdate(double X, double_indexed *Y, int fold) {
  dmdupdate(X, Y, 1, Y + fold, 1, fold);
}

void zmdupdate(double X, double* repY, int increpY, double* carY, int inccarY, int fold) {
  dmdupdate(X, repY, 2 * increpY, carY, 2 * inccarY, fold);
  dmdupdate(X, repY + 1, 2 * increpY, carY + 1, 2 * inccarY, fold);
}

void zidupdate(double X, double_complex_indexed *Y, int fold) {
  zmdupdate(X, Y, 1, Y + 2 * fold, 1, fold);
}

void zmzupdate(void *X, double* repY, int increpY, double* carY, int inccarY, int fold) {
  dmdupdate(((double*)X)[0], repY, 2 * increpY, carY, 2 * inccarY, fold);
  dmdupdate(((double*)X)[1], repY + 1, 2 * increpY, carY + 1, 2 * inccarY, fold);
}

void zizupdate(void *X, double_complex_indexed *Y, int fold) {
  zmzupdate(X, Y, 1, Y + 2 * fold, 1, fold);
}
