#include <math.h>

#include <indexed.h>

#include "../common/common.h"

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
  int j;
  const float *bins;

  if (isnan(manY[0]) || isinf(manY[0]))
    return;

  int X_index = sindex(X);
  int shift = smindex(manY) - X_index;
  if(shift > 0){
    for(i = fold - 1; i >= shift; i--){
      manY[i * incmanY] = manY[(i - shift) * incmanY];
      carY[i * inccarY] = carY[(i - shift) * inccarY];
    }
    bins = smbins(X_index);
    for(j = 0; j < i + 1; j++){
      manY[j * incmanY] = bins[j];
      carY[j * inccarY] = 0.0;
    }
  }
}
