#include <indexed.h>

#include <math.h>

#include "../common/common.h"

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
  int i = 0;
  double Y = 0.0;
  int X_index;
  const float *bins;

  if (ISNANINFF(manX[0])){
    return manX[0];
  }

  if (manX[0] == 0.0) {
    return 0.0;
  }

  //Note that the following order of summation is in order of decreasing
  //exponent. The following code is specific to SIWIDTH=13, FLT_MANT_DIG=24, and
  //the number of carries equal to 1.
  X_index = smindex(manX);
  bins = smbins(X_index);
  if(X_index == 0){
    Y += (double)carX[0] * (double)(bins[0]/6.0) * (double)SMEXPANSION;
    if(fold > 1){
      Y += (double)carX[inccarX] * (double)(bins[1]/6.0);
    }
    Y += (double)(manX[0] - bins[0]) * (double)SMEXPANSION;
    if(fold == 1){
      return (float)Y;
    }
    i = 2;
  }else{
    Y += (double)carX[0] * (double)(bins[0]/6.0);
    i = 1;
  }
  for(; i < fold; i++){
    Y += (double)carX[i * inccarX] * (double)(bins[i]/6.0);
    Y += (double)(manX[(i - 1) * incmanX] - bins[i - 1]);
  }
  Y += (double)(manX[(fold - 1) * incmanX] - bins[fold - 1]);

  return (float)Y;
}
