/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

void smsupdate(float X, float* repY, int increpY, float* carY, int inccarY, int fold) {
  if (X == 0 || isnan(repY[0]) || isinf(repY[0]))
    return;

/*
  if (repY[0] == 0.0) {
    smbound(sindex(X), repY, increpY, fold);
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
    smbound(X_index, repY, increpY, MIN(d, fold));
    for(int i = 0; i < d && i < fold; i++){
      carY[i * inccarY] = 0.0;
    }
  }
}

void sisupdate(float X, float_indexed *Y, int fold) {
  smsupdate(X, Y, 1, Y + fold, 1, fold);
}

void cmsupdate(float X, float* repY, int increpY, float* carY, int inccarY, int fold) {
  smsupdate(X, repY, 2 * increpY, carY, 2 * inccarY, fold);
  smsupdate(X, repY + 1, 2 * increpY, carY + 1, 2 * inccarY, fold);
}

void cisupdate(float X, float_complex_indexed *Y, int fold) {
  cmsupdate(X, Y, 1, Y + 2 * fold, 1, fold);
}

void cmcupdate(void *X, float* repY, int increpY, float* carY, int inccarY, int fold) {
  smsupdate(((float*)X)[0], repY, 2 * increpY, carY, 2 * inccarY, fold);
  smsupdate(((float*)X)[1], repY + 1, 2 * increpY, carY + 1, 2 * inccarY, fold);
}

void cicupdate(void *X, float_complex_indexed *Y, int fold) {
  cmcupdate(X, Y, 1, Y + 2 * fold, 1, fold);
}
