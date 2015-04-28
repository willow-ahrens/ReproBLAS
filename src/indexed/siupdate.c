/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

/**
 * @brief Update manually specified indexed single precision with single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value less than X
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void smsupdate(const int fold, float X, float* repY, int increpY, float* carY, int inccarY) {
  if (X == 0 || isnan(repY[0]) || isinf(repY[0]))
    return;

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

/**
 * @brief Update indexed single precision with single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value less than X
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void sisupdate(const int fold, float X, float_indexed *Y) {
  smsupdate(fold, X, Y, 1, Y + fold, 1);
}

/**
 * @brief Update manually specified indexed complex single precision with single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value less than X
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmsupdate(const int fold, float X, float* repY, int increpY, float* carY, int inccarY) {
  smsupdate(fold, X, repY, 2 * increpY, carY, 2 * inccarY);
  smsupdate(fold, X, repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}

/**
 * @brief Update indexed complex single precision with single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value less than X
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cisupdate(const int fold, float X, float_complex_indexed *Y) {
  cmsupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}

/**
 * @brief Update manually specified indexed complex single precision with complex single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value of real and imaginary components less than absolute value of real and imaginary components of X respectively.
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmcupdate(const int fold, void *X, float* repY, int increpY, float* carY, int inccarY) {
  smsupdate(fold, ((float*)X)[0], repY, 2 * increpY, carY, 2 * inccarY);
  smsupdate(fold, ((float*)X)[1], repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}

/**
 * @brief Update indexed complex single precision with complex single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value of real and imaginary components less than absolute value of real and imaginary components of X respectively.
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cicupdate(const int fold, void *X, float_complex_indexed *Y) {
  cmcupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}
