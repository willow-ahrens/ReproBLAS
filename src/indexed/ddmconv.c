#include <math.h>

#include <indexed.h>

#include "../../config.h"

#include "../common/common.h"

void ddpd(double* a, double b) {
  (void)ddpd;
  double q;
  double s1, s2, t1, t2;

  /* Add two hi words. */
  s1 = a[0] + b;
  q = s1 - a[0];
  s2 = ((b - q) + (a[0] - (s1 - q)));

  t1 = a[1] + s2;
  q = t1 - a[1];
  t2 = ((s2 - q) + (a[1] - (t1 - q)));

  s2 = t1;

  /* Renormalize (s1, s2)  to  (t1, s2) */
  t1 = s1 + s2;
  t2 += s2 - (t1 - s1);

  /* Renormalize (t1, t2)  */
  a[0] = t1 + t2;
  a[1] = t2 - (a[0] - t1);
}

/* Purpose
 * =======
 *
 * This subroutine computes ddc = dda + ddb.
 *
 * Taken from D. H. Bailey's ddfun90.f.
 *
 */
/*
void ddadd(double dda_l, double dda_t, double ddb_l, double ddb_t,
           double *ddc_l, double *ddc_t)
{
  double e, t1, t2;

  // Compute dda + ddb using Knuth's trick.
  t1 = dda_l + ddb_l;
  e = t1 - dda_l;
  t2 = ((ddb_l - e) + (dda_l - (t1 - e))) + dda_t + ddb_t;

  // The result is t1 + t2, after normalization.
  *ddc_l = t1 + t2;
  *ddc_t = t2 - (*ddc_l - t1);

}
*/

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

  X_index = dmindex(manX);
  bins = dmbins(X_index);

  //Note that the following order of summation is in order of decreasing
  //exponent. The following code is specific to DIWIDTH=13, DBL_MANT_DIG=24, and
  //the number of carries equal to 1.
  #if (LDBL_MANT_DIG - DBL_MANT_DIG) > 15 || 1 << (LDBL_MANT_DIG - DBL_MANT_DIG) >= 2 * MAX_FOLD
    #if LDBL_MAX_EXP > DBL_MAX_EXP + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2)
      long double Y = 0.0;
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
      return (double)Y;

    #else
      long double Y = 0.0;
      double scale_down = ldexp(0.5, -1 * (DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2 + 1)) + 1);
      double scale_up = ldexp(0.5, DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2 + 1) + 1);
      X_index = dmindex(manX);
      bins = dmbins(X_index);
      int scaled = MIN(X_index + fold, (DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2))/DIWIDTH);
      if(X_index <= (DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2))/DIWIDTH){
        if(X_index == 0){
          Y += carX[0] * ((bins[0]/6.0) * scale_down * DMEXPANSION);
          if(fold > 1){
            Y += carX[inccarX] * ((bins[1]/6.0) * scale_down);
          }
          Y += (manX[0] - bins[0]) * scale_down * DMEXPANSION;
          i = 2;
        }else{
          Y += carX[0] * ((bins[0]/6.0) * scale_down);
          i = 1;
        }
        for(; i < scaled && i < fold; i++){
          Y += carX[i * inccarX] * ((bins[i]/6.0) * scale_down);
          Y += (manX[(i - 1) * incmanX] - bins[i - 1]) * scale_down;
        }
        if(i == fold){
          Y += (manX[(fold - 1) * incmanX] - bins[fold - 1]) * scale_down;
          return Y * scale_up;
        }
        Y *= scale_up;
        for(; i < fold; i++){
          Y += carX[i * inccarX] * (bins[i]/6.0);
          Y += manX[(i - 1) * incmanX] - bins[i - 1];
        }
        Y += manX[(fold - 1) * incmanX] - bins[fold - 1];
        return (double)Y;
      }else{
        Y += carX[0] * (bins[0]/6.0);
        for(i = 1; i < fold; i++){
          Y += carX[i * inccarX] * (bins[i]/6.0);
          Y += (manX[(i - 1) * incmanX] - bins[i - 1]);
        }
        Y += (manX[(fold - 1) * incmanX] - bins[fold - 1]);
        return (double)Y;
      }

    #endif
  #else
    double Y[2] = {0.0, 0.0};
    double scale_down = ldexp(0.5, -1 * (DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2 + 1)) + 1);
    double scale_up = ldexp(0.5, DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2 + 1) + 1);
    X_index = dmindex(manX);
    bins = dmbins(X_index);
    int scaled = MIN(X_index + fold, (DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2))/DIWIDTH);
    if(X_index <= (DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2))/DIWIDTH){
      if(X_index == 0){
        ddpd(Y, carX[0] * ((bins[0]/6.0) * scale_down * DMEXPANSION));
        if(fold > 1){
          ddpd(Y, carX[inccarX] * ((bins[1]/6.0) * scale_down));
        }
        ddpd(Y, (manX[0] - bins[0]) * scale_down * DMEXPANSION);
        i = 2;
      }else{
        ddpd(Y, carX[0] * ((bins[0]/6.0) * scale_down));
        i = 1;
      }
      for(; i < scaled && i < fold; i++){
        ddpd(Y, carX[i * inccarX] * ((bins[i]/6.0) * scale_down));
        ddpd(Y, (manX[(i - 1) * incmanX] - bins[i - 1]) * scale_down);
      }
      if(i == fold){
        ddpd(Y, (manX[(fold - 1) * incmanX] - bins[fold - 1]) * scale_down);
        return (Y[0] + Y[1]) * scale_up;
      }
      Y[0] *= scale_up;
      Y[1] *= scale_up;
      for(; i < fold; i++){
        ddpd(Y, carX[i * inccarX] * (bins[i]/6.0));
        ddpd(Y, manX[(i - 1) * incmanX] - bins[i - 1]);
      }
      ddpd(Y, manX[(fold - 1) * incmanX] - bins[fold - 1]);
      return Y[0] + Y[1];
    }else{
      ddpd(Y, carX[0] * (bins[0]/6.0));
      for(i = 1; i < fold; i++){
        ddpd(Y, carX[i * inccarX] * (bins[i]/6.0));
        ddpd(Y, (manX[(i - 1) * incmanX] - bins[i - 1]));
      }
      ddpd(Y, (manX[(fold - 1) * incmanX] - bins[fold - 1]));
      return Y[0] + Y[1];
    }
  #endif
}
