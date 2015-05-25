#include <math.h>

#include "indexed.h"
#include "../Common/Common.h"

/**
 * @internal
 * @brief Update manually specified indexed single precision with single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value less than X
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   5 May 2015
 */
void smsupdate(const int fold, const float X, float* manY, const int incmanY, float* carY, const int inccarY) {
  int i;

  if (X == 0 || isnan(manY[0]) || isinf(manY[0]))
    return;

  int X_index = sindex(X);
  int shift = smindex(manY) - X_index;
  if(shift > 0){
    for(i = fold - 1; i >= shift; i--){
      manY[i * incmanY] = manY[(i - shift) * incmanY];
      carY[i * inccarY] = carY[(i - shift) * inccarY];
    }
    smbin(MIN(shift, fold), X_index, manY, incmanY, carY, inccarY);
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
void sisupdate(const int fold, const float X, float_indexed *Y) {
  smsupdate(fold, X, Y, 1, Y + fold, 1);
}

/**
 * @internal
 * @brief Update manually specified indexed complex single precision with single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value less than X
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmsupdate(const int fold, const float X, float* manY, const int incmanY, float* carY, const int inccarY) {
  smsupdate(fold, X, manY, 2 * incmanY, carY, 2 * inccarY);
  smsupdate(fold, X, manY + 1, 2 * incmanY, carY + 1, 2 * inccarY);
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
void cisupdate(const int fold, const float X, float_complex_indexed *Y) {
  cmsupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}

/**
 * @internal
 * @brief Update manually specified indexed complex single precision with complex single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value of real and imaginary components less than absolute value of real and imaginary components of X respectively.
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmcupdate(const int fold, const void *X, float* manY, const int incmanY, float* carY, const int inccarY) {
  smsupdate(fold, ((float*)X)[0], manY, 2 * incmanY, carY, 2 * inccarY);
  smsupdate(fold, ((float*)X)[1], manY + 1, 2 * incmanY, carY + 1, 2 * inccarY);
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
void cicupdate(const int fold, const void *X, float_complex_indexed *Y) {
  cmcupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}
