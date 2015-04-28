/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

void smsupdate(const int fold, float X, float* repY, int increpY, float* carY, int inccarY) {
  if (X == 0 || isnan(repY[0]) || isinf(repY[0]))
    return;

/*
  if (repY[0] == 0.0) {
    smbound(fold, sindex(X), repY, increpY);
    for (int i = 0; i < fold; i++) {
      carY[i * inccarY] = 0.0;
    }
    return;
  }
*/

  int X_index = sindex(X);
  int d = smindex(repY) - X_index;
  if(d > 0){
    for(int i = fold - 1; i >= d; i--){
      repY[i * increpY] = repY[(i - d) * increpY];
      carY[i * inccarY] = carY[(i - d) * inccarY];
    }
    smbound(MIN(d, fold), X_index, repY, increpY);
    for(int i = 0; i < d && i < fold; i++){
      carY[i * inccarY] = 0.0;
    }
  }
}

void sisupdate(const int fold, float X, float_indexed *Y) {
  smsupdate(fold, X, Y, 1, Y + fold, 1);
}

void cmsupdate(const int fold, float X, float* repY, int increpY, float* carY, int inccarY) {
  smsupdate(fold, X, repY, 2 * increpY, carY, 2 * inccarY);
  smsupdate(fold, X, repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}

void cisupdate(const int fold, float X, float_complex_indexed *Y) {
  cmsupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}

void cmcupdate(const int fold, void *X, float* repY, int increpY, float* carY, int inccarY) {
  smsupdate(fold, ((float*)X)[0], repY, 2 * increpY, carY, 2 * inccarY);
  smsupdate(fold, ((float*)X)[1], repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}

void cicupdate(const int fold, void *X, float_complex_indexed *Y) {
  cmcupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}
