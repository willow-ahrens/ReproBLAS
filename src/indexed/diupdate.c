/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

/**
 * @internal
 * @brief Update manually specified indexed double precision with double precision (X -> Y)
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
void dmdupdate(const int fold, const double X, double* repY, const int increpY, double* carY, const int inccarY) {
  if (X == 0 || isnan(repY[0]) || isinf(repY[0]))
    return;

  int X_index = dindex(X);
  int shift = dmindex(repY) - X_index;
  if(shift > 0){
    for(int i = fold - 1; i >= shift; i--){
      repY[i * increpY] = repY[(i - shift) * increpY];
      carY[i * inccarY] = carY[(i - shift) * inccarY];
    }
    dmbound(MIN(shift, fold), X_index, repY, increpY, carY, inccarY);
  }
}

/**
 * @brief Update indexed double precision with double precision (X -> Y)
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
void didupdate(const int fold, const double X, double_indexed *Y) {
  dmdupdate(fold, X, Y, 1, Y + fold, 1);
}

/**
 * @internal
 * @brief Update manually specified indexed complex double precision with double precision (X -> Y)
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
void zmdupdate(const int fold, const double X, double* repY, const int increpY, double* carY, const int inccarY) {
  dmdupdate(fold, X, repY, 2 * increpY, carY, 2 * inccarY);
  dmdupdate(fold, X, repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}

/**
 * @brief Update indexed complex double precision with double precision (X -> Y)
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
void zidupdate(const int fold, const double X, double_complex_indexed *Y) {
  zmdupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}

/**
 * @internal
 * @brief Update manually specified indexed complex double precision with complex double precision (X -> Y)
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
void zmzupdate(const int fold, const void *X, double* repY, const int increpY, double* carY, const int inccarY) {
  dmdupdate(fold, ((double*)X)[0], repY, 2 * increpY, carY, 2 * inccarY);
  dmdupdate(fold, ((double*)X)[1], repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}

/**
 * @brief Update indexed complex double precision with complex double precision (X -> Y)
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
void zizupdate(const int fold, const void *X, double_complex_indexed *Y) {
  zmzupdate(fold, X, Y, 1, Y + 2 * fold, 1);
}
