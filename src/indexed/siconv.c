#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

/**
 * @internal
 * @brief Convert single precision to manually specified indexed single precision (X -> Y)
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
void smsconv(const int fold, const float X, float* manY, const int incmanY, float* carY, const int inccarY) {
  smsetzero(fold, manY, incmanY, carY, inccarY);
  smsadd(fold, X, manY, incmanY, carY, inccarY);
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
void sisconv(const int fold, const float X, float_indexed *Y) {
  smsconv(fold, X, Y, 1, Y + fold, 1);
}

/**
 * @internal
 * @brief Convert complex single precision to manually specified indexed complex single precision (X -> Y)
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
void cmcconv(const int fold, const void *X, float *manY, const int incmanY, float *carY, const int inccarY) {
  smsconv(fold, ((float*)X)[0], manY, incmanY * 2, carY, inccarY * 2);
  smsconv(fold, ((float*)X)[1], manY + 1, incmanY * 2, carY + 1, inccarY * 2);
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
void cicconv(const int fold, const void *X, float_complex_indexed *Y) {
  cmcconv(fold, X, Y, 1, Y + 2 * fold, 1);
}

/**
 * @internal
 * @brief Convert manually specified indexed single precision to single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @return scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
float ssmconv(const int fold, const float* manX, const int incmanX, const float* carX, const int inccarX) {
  int i;
  float Y = 0.0;

  if (isinf(manX[0]) || isnan(manX[0]))
    return manX[0];

  if (manX[0] == 0.0) {
    return 0.0;
  }

  // TODO: SCALING TO AVOID OVERFLOW

  for (i = 0; i < fold; i++) {
    Y += (manX[i * incmanX] + (carX[i * incmanX] - 6) * ufpf(manX[i * incmanX]) * 0.25);
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
float ssiconv(const int fold, const float_indexed *X) {
  return ssmconv(fold, X, 1, X + fold, 1);
}

/**
 * @internal
 * @brief Convert manually specified indexed complex single precision to complex single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param conv scalar return
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void ccmconv_sub(const int fold, const float *manX, const int incmanX, const float *carX, const int inccarX, void *conv) {
  ((float*)conv)[0] = ssmconv(fold, manX, incmanX * 2, carX, inccarX + 1);
  ((float*)conv)[1] = ssmconv(fold, manX + 1, incmanX * 2, carX + 1, inccarX + 1);
}

/**
 * @brief Convert indexed complex single precision to complex single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param conv scalar return
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cciconv_sub(const int fold, const float_complex_indexed *X, void *conv) {
  ccmconv_sub(fold, X, 1, X + 2 * fold, 1, conv);
}
