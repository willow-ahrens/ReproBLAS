#include <math.h>

#include <indexed.h>

#include "../../config.h"

/**
 * @internal
 * @brief Get index of manually specified indexed double precision
 *
 * The index of an indexed type is the bin that it corresponds to. Higher indicies correspond to smaller bins.
 *
 * @param manX X's mantissa vector
 * @return X's index
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   19 May 2015
 */
int dmindex(const double *manX){
  int exp;

  if(manX[0] == 0.0){
    return (DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH + MAX_FOLD;
  }else{
    frexp(manX[0], &exp);
    if(exp == DBL_MAX_EXP){
      return 0;
    }
    return (DBL_MAX_EXP + DBL_MANT_DIG - DIWIDTH + 1 - exp)/DIWIDTH;
  }
  //return (DBL_MAX_EXP + DBL_MANT_DIG - DIWIDTH + 1 - EXP(manX[0]))/DIWIDTH;
}
