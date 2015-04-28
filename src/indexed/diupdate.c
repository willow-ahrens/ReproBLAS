/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

void dmdupdate(const int fold, double X, double* repY, int increpY, double* carY, int inccarY) {
  if (X == 0 || isnan(repY[0]) || isinf(repY[0]))
    return;

/*
  if (repY[0] == 0.0) {
    dmbound(fold, dindex(X), repY, increpY);
    for (int i = 0; i < fold; i++) {
      carY[i * inccarY] = 0.0;
    }
    return;
  }
*/

  int X_index = dindex(X);
  int d = dmindex(repY) - X_index;
  if(d > 0){
    for(int i = fold - 1; i >= d; i--){
      repY[i * increpY] = repY[(i - d) * increpY];
      carY[i * inccarY] = carY[(i - d) * inccarY];
    }
    dmbound(MIN(d, fold), X_index, repY, increpY);
    for(int i = 0; i < d && i < fold; i++){
      carY[i * inccarY] = 0.0;
    }
  }
}

void didupdate(const int fold, double X, double_indexed *Y) {
  dmdupdate(fold, X, Y, 1, Y + fold, 1);
}

void zmdupdate(const int fold, double X, double* repY, int increpY, double* carY, int inccarY) {
  dmdupdate(fold, X, repY, 2 * increpY, carY, 2 * inccarY);
  dmdupdate(fold, X, repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}

void zidupdate(const int fold, double X, double_complex_indexed *Y) {
  zmdupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}

void zmzupdate(const int fold, void *X, double* repY, int increpY, double* carY, int inccarY) {
  dmdupdate(fold, ((double*)X)[0], repY, 2 * increpY, carY, 2 * inccarY);
  dmdupdate(fold, ((double*)X)[1], repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}

void zizupdate(const int fold, void *X, double_complex_indexed *Y) {
  zmzupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}
