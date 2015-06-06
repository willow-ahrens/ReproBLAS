#include <math.h>

#include <indexed.h>

#include <../../config.h>

/**
 * @internal
 * @brief Convert manually specified indexed double precision to double precision (X -> Y)
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
double ddmconv(const int fold, const double* manX, const int incmanX, const double* carX, const int inccarX) {
  int i = 0;
  int X_index;
  long double Y = 0.0;
  const double *bins;

  if (isinf(manX[0]) || isnan(manX[0]))
    return manX[0];

  if (manX[0] == 0.0) {
    return 0.0;
  }

  /*
  //Naive
  double M;
  if(dmindex0(manX)){
    M = ufp(manX[i * incmanX]);
    Y += (carX[i * inccarX] * 0.25 * M + (manX[i * incmanX] - 1.5 * M)) * DMEXPANSION;
    i = 1;
  }

  for (; i < fold; i++) {
    M = ufp(manX[i * incmanX]);
    Y += carX[i * inccarX] * 0.25 * M + (manX[i * incmanX] - 1.5 * M);
  }
  */

  //TODO Double check that max carry value is 2^53

  //Note that the following order of summation is in order of decreasing
  //exponent. The following code is specific to DIWIDTH=13, DBL_MANT_DIG=24, and
  //the number of carries equal to 1.
  #if (LDBL_MANT_DIG - DBL_MANT_DIG) > 15 || 1 << (LDBL_MANT_DIG - DBL_MANT_DIG) > (MAX_FLOW + 1) * MAX_FOLD
    #if LDBL_MAX_EXP > DBL_MAX_EXP + (DBL_MANT_DIG * MAX_FLOW) + (DBL_MANT_DIG - DIWIDTH - 2)
      X_index = dmindex(manX);
      bins = dmbins(X_index);
      if(X_index == 0){
        Y += (long double)carX[0] * (long double)(bins[0]/6.0) * (long double)DMEXPANSION;
        if(fold > 1){
          Y += (long double)carX[inccarX] * (long double)(bins[1]/6.0);
        }
        Y += (long double)(manX[0] - bins[0]) * (long double)DMEXPANSION;
        i = 2;
      }else{
        Y += (long double)carX[0] * (long double)(bins[0]/6.0);
        i = 1;
      }
      for(; i < fold; i++){
        Y += (long double)carX[i * inccarX] * (long double)(bins[i]/6.0);
        Y += (long double)(manX[(i - 1) * incmanX] - bins[i - 1]);
      }
      Y += (long double)(manX[(fold - 1) * incmanX] - bins[fold - 1]);
    #else
      not actually code
    #endif
  #else
    not actually code
  #endif

  return (double)Y;
}
