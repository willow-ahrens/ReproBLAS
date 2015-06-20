#include <math.h>

#include <indexed.h>

#include "../common/common.h"

#include "../../config.h"

#define DDPD(T0, T1, T2, T3, T4, T5, Yl, Yt, X) \
  T0 = (X); \
  T1 = Yl + T0; \
  T2 = T1 - Yl; \
  T3 = ((T0 - T2) + (Yl - (T1 - T2))); \
  T4 = Yt + T3; \
  T2 = T4 - Yt; \
  T5 = ((T3 - T2) + (Yt - (T4 - T2))); \
  T3 = T4; \
  T4 = T1 + T3; \
  T5 += T3 - (T4 - T1); \
  Yl = T4 + T5; \
  Yt = T5 - (Yl - T4);

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

  if (ISNANINF(manX[0])){
    return manX[0];
  }

  if (manX[0] == 0.0) {
    return 0.0;
  }

  #if 0
    //Naive
    double Y = 0.0;
    X_index = dmindex(manX);
    bins = dmbins(X_index);
    if(X_index == 0){
      Y += (carX[i * inccarX] * (bins[0]/6.0) + (manX[i * incmanX] - bins[0])) * DMEXPANSION;
      i = 1;
    }

    for (; i < fold; i++) {
      Y += carX[i * inccarX] * (bins[i]/6.0) + (manX[i * incmanX] - bins[i]);
    }
    return Y;
  #else
    /*
    */
    //Note that the following order of summation is in order of decreasing
    //exponent. The following code is specific to DIWIDTH=13, DBL_MANT_DIG=24, and
    //the number of carries equal to 1.
    X_index = dmindex(manX);
    bins = dmbins(X_index);
    #if 0 && ((LDBL_MANT_DIG - DBL_MANT_DIG) > 15 || 1 << (LDBL_MANT_DIG - DBL_MANT_DIG) >= 2 * MAX_FOLD)
      #if 0 && LDBL_MAX_EXP > DBL_MAX_EXP + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2)
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
        double Y = 0.0;
        double scale_down;
        double scale_up;
        int scaled;
        X_index = dmindex(manX);
        bins = dmbins(X_index);
        if(X_index <= (DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2))/DIWIDTH){
          scale_down = ldexp(0.5, 1 - (DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2 + 1)));
          scale_up = ldexp(0.5, 1 + DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2 + 1));
          scaled = MIN(X_index + fold, (DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2))/DIWIDTH);
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
          if(isinf(Y)){
            return Y;
          }
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
      double Yl = 0.0;
      double Yt = 0.0;
      double t0;
      double t1;
      double t2;
      double t3;
      double t4;
      double t5;
      double scale_down;
      double scale_up;
      int scaled;
      X_index = dmindex(manX);
      bins = dmbins(X_index);
      if(X_index <= (DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2))/DIWIDTH){
        scale_down = ldexp(0.5, 1 - (DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2 + 1)));
        scale_up = ldexp(0.5, 1 + DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2 + 1));
        scaled = MIN(X_index + fold, (DBL_MANT_DIG + (DBL_MANT_DIG - 1) + (DBL_MANT_DIG - DIWIDTH - 2))/DIWIDTH);
        if(X_index == 0){
          DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, carX[0] * ((bins[0]/6.0) * scale_down * DMEXPANSION));
          DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, carX[inccarX] * ((bins[1]/6.0) * scale_down));
          DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, (manX[0] - bins[0]) * scale_down * DMEXPANSION);
          i = 2;
        }else{
          DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, carX[0] * ((bins[0]/6.0) * scale_down));
          i = 1;
        }
        for(; i < scaled && i < fold; i++){
          DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, carX[i * inccarX] * ((bins[i]/6.0) * scale_down));
          DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, (manX[(i - 1) * incmanX] - bins[i - 1]) * scale_down);
        }
        if(i == fold){
          DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, (manX[(fold - 1) * incmanX] - bins[fold - 1]) * scale_down);
          return (Yl + Yt) * scale_up;
        }
        if(isinf((Yl + Yt) * scale_up)){
          return (Yl + Yt) * scale_up;
        }
        Yl *= scale_up;
        Yt *= scale_up;
        for(; i < fold; i++){
          DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, carX[i * inccarX] * (bins[i]/6.0));
          DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, manX[(i - 1) * incmanX] - bins[i - 1]);
        }
        DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, manX[(fold - 1) * incmanX] - bins[fold - 1]);
        return Yl + Yt;
      }else{
        DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, carX[0] * (bins[0]/6.0));
        for(i = 1; i < fold; i++){
          DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, carX[i * inccarX] * (bins[i]/6.0));
          DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, (manX[(i - 1) * incmanX] - bins[i - 1]));
        }
        DDPD(t0, t1, t2, t3, t4, t5, Yl, Yt, (manX[(fold - 1) * incmanX] - bins[fold - 1]));
        return Yl + Yt;
      }
    #endif
  #endif
}
