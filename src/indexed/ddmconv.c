#include <math.h>

#include <indexed.h>

#include "../common/common.h"

#include "../../config.h"

/**
 * @internal
 * @brief Convert manually specified indexed double precision to double precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @return scalar Y
 *
 * @author Peter Ahrens
 * @date   31 Jul 2015
 */
double ddmconv(const int fold, const double* priX, const int incpriX, const double* carX, const int inccarX) {
  int i = 0;
  int X_index;
  const double *bins;

  if (ISNANINF(priX[0])){
    return priX[0];
  }

  if (priX[0] == 0.0) {
    return 0.0;
  }

  double Y = 0.0;
  double scale_down;
  double scale_up;
  int scaled;
  X_index = dmindex(priX);
  bins = dmbins(X_index);
  if(X_index <= (3 * DBL_MANT_DIG - 1)/DIWIDTH){
    scale_down = ldexp(0.5, 1 - (2 * DBL_MANT_DIG - DIWIDTH - 2));
    scale_up = ldexp(0.5, 1 + (2 * DBL_MANT_DIG - DIWIDTH - 2));
    scaled = MAX(MIN(fold, (3 * DBL_MANT_DIG - 1)/DIWIDTH - X_index), 0);
    if(X_index == 0){
      Y += carX[0] * ((bins[0]/6.0) * scale_down * DMEXPANSION);
      Y += carX[inccarX] * ((bins[1]/6.0) * scale_down);
      Y += (priX[0] - bins[0]) * scale_down * DMEXPANSION;
      i = 2;
    }else{
      Y += carX[0] * ((bins[0]/6.0) * scale_down);
      i = 1;
    }
    for(; i < scaled; i++){
      Y += carX[i * inccarX] * ((bins[i]/6.0) * scale_down);
      Y += (priX[(i - 1) * incpriX] - bins[i - 1]) * scale_down;
    }
    if(i == fold){
      Y += (priX[(fold - 1) * incpriX] - bins[fold - 1]) * scale_down;
      return Y * scale_up;
    }
    if(isinf(Y * scale_up)){
      return Y * scale_up;
    }
    Y *= scale_up;
    for(; i < fold; i++){
      Y += carX[i * inccarX] * (bins[i]/6.0);
      Y += priX[(i - 1) * incpriX] - bins[i - 1];
    }
    Y += priX[(fold - 1) * incpriX] - bins[fold - 1];
  }else{
    Y += carX[0] * (bins[0]/6.0);
    for(i = 1; i < fold; i++){
      Y += carX[i * inccarX] * (bins[i]/6.0);
      Y += (priX[(i - 1) * incpriX] - bins[i - 1]);
    }
    Y += (priX[(fold - 1) * incpriX] - bins[fold - 1]);
  }
  return Y;
}
