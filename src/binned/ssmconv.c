#include <binned.h>

#include <math.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Convert manually specified binned single precision to single precision (X -> Y)
 *
 * @param fold the fold of the binned types
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @return scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
float binned_ssmconv(const int fold, const float* priX, const int incpriX, const float* carX, const int inccarX) {
  int i = 0;
  double Y = 0.0;
  int X_index;
  const float *bins;

  if (ISNANINFF(priX[0])){
    return priX[0];
  }

  if (priX[0] == 0.0) {
    return 0.0;
  }

  //Note that the following order of summation is in order of decreasing
  //exponent. The following code is specific to SBWIDTH=13, FLT_MANT_DIG=24, and
  //the number of carries equal to 1.
  X_index = binned_smindex(priX);
  bins = binned_smbins(X_index);
  if(X_index == 0){
    Y += (double)carX[0] * (double)(bins[0]/6.0) * (double)binned_SMEXPANSION;
    Y += (double)carX[inccarX] * (double)(bins[1]/6.0);
    Y += (double)(priX[0] - bins[0]) * (double)binned_SMEXPANSION;
    i = 2;
  }else{
    Y += (double)carX[0] * (double)(bins[0]/6.0);
    i = 1;
  }
  for(; i < fold; i++){
    Y += (double)carX[i * inccarX] * (double)(bins[i]/6.0);
    Y += (double)(priX[(i - 1) * incpriX] - bins[i - 1]);
  }
  Y += (double)(priX[(fold - 1) * incpriX] - bins[fold - 1]);

  return (float)Y;
}
