#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

/**
 * @brief Convert single precision to manually specified indexed single precision (X -> Y)
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
void smsconv(const int fold, float X, float* repY, int increpY, float* carY, int inccarY) {
  int i;
  float q;
  float s;
  if (X == 0.0) {
    for (i = 0; i < fold; i++) {
      repY[i*increpY] = 0.0;
      carY[i*inccarY] = 0.0;
    }
    return;
  }
  smbound(fold, sindex(fabs(X)), repY, increpY);
  for (i = 0; i < fold; i++, repY += increpY, carY += inccarY) {
    carY[0] = 0.0;

    s = repY[0];
    q = s + X;
    repY[0] = s;
    q -= s;
    X -= q;
  }
}

/**
 * @brief Convert single precision to indexed single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void sisconv(const int fold, float X, float_indexed *Y) {
  smsconv(fold, X, Y, 1, Y + fold, 1);
}

/**
 * @brief Convert complex single precision to manually specified indexed complex single precision (X -> Y)
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
void cmcconv(const int fold, void *X, float *repY, int increpY, float *carY, int inccarY) {
  smsconv(fold, ((float*)X)[0], repY, increpY * 2, carY, inccarY * 2);
  smsconv(fold, ((float*)X)[1], repY + 1, increpY * 2, carY + 1, inccarY * 2);
}

/**
 * @brief Convert complex single precision to indexed complex single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cicconv(const int fold, void *X, float_complex_indexed *Y) {
  cmcconv(fold, X, Y, 1, Y + 2 * fold, 1);
}

/**
 * @brief Convert manually specified indexed single precision to single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param repX X's rep vector
 * @param increpX stride within X's rep vector (use every increpX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @return scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
float ssmconv(const int fold, float* repX, int increpX, float* carX, int inccarX) {
  int i;
  float Y = 0.0;

  if (isinf(repX[0]) || isnan(repX[0]))
    return repX[0];

  if (repX[0] == 0.0) {
    return 0.0;
  }

  // TODO: SCALING TO AVOID OVERFLOW

  for (i = 0; i < fold; i++, repX += increpX, carX += inccarX) {
    Y += (repX[0] + (carX[0] - 6) * ufp(repX[0]) * 0.25);
  }

  return Y;
}

/**
 * @brief Convert indexed single precision to single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @return scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
float ssiconv(const int fold, float_indexed *X) {
  return ssmconv(fold, X, 1, X + fold, 1);
}

/**
 * @brief Convert manually specified indexed complex single precision to complex single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param repX X's rep vector
 * @param increpX stride within X's rep vector (use every increpX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param Y scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void ccmconv_sub(const int fold, float *repX, int increpX, float *carX, int inccarX, void *Y) {
  ((float*)Y)[0] = ssmconv(fold, repX, increpX * 2, carX, inccarX + 1);
  ((float*)Y)[1] = ssmconv(fold, repX + 1, increpX * 2, carX + 1, inccarX + 1);
}

/**
 * @brief Convert indexed complex single precision to complex single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param Y scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cciconv_sub(const int fold, float_complex_indexed *X, void *Y) {
  ccmconv_sub(fold, X, 1, X + 2 * fold, 1, Y);
}
