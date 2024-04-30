#include <math.h>

#include <binned.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Update manually specified binned single precision with single precision (X -> Y)
 *
 * This method updates Y to an index suitable for adding numbers with absolute value less than X
 *
 * @param fold the fold of the binned types
 * @param X scalar X
 * @param priY Y's primary vector
 * @param incpriY stride within Y's primary vector (use every incpriY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   5 May 2015
 */
void binned_smsupdate(const int fold, const float X, float* priY, const int incpriY, float* carY, const int inccarY) {
  int i;
  int j;
  int X_index;
  int shift;
  const float *bins;

  if (ISNANINFF(priY[0])){
    return;
  }

  X_index = binned_sindex(X);
  if(priY[0] == 0.0){
    bins = binned_smbins(X_index);
    for(i = 0; i < fold; i++){
      priY[i * incpriY] = bins[i];
      carY[i * inccarY] = 0.0;
    }
  }else{
    shift = binned_smindex(priY) - X_index;
    if(shift > 0){
      for(i = fold - 1; i >= shift; i--){
        priY[i * incpriY] = priY[(i - shift) * incpriY];
        carY[i * inccarY] = carY[(i - shift) * inccarY];
      }
      bins = binned_smbins(X_index);
      for(j = 0; j < i + 1; j++){
        priY[j * incpriY] = bins[j];
        carY[j * inccarY] = 0.0;
      }
    }
  }
}
