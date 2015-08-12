#include <math.h>

#include <idxd.h>

#include "../common/common.h"

#include "../../config.h"

/**
 * @internal
 * @brief Get index of manually specified indexed double precision
 *
 * The index of an indexed type is the bin that it corresponds to. Higher indicies correspond to smaller bins.
 *
 * @param priX X's primary vector
 * @return X's index
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   19 May 2015
 */
int idxd_dmindex(const double *priX){
  /*
  int exp;

  if(priX[0] == 0.0){
    return (DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH + DIMAXFOLD;
  }else{
    frexp(priX[0], &exp);
    if(exp == DBL_MAX_EXP){
      return 0;
    }
    return (DBL_MAX_EXP + DBL_MANT_DIG - DIWIDTH + 1 - exp)/DIWIDTH;
  }
  */
  /*
  if(priX[0] == 0.0){
    return (DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH + DIMAXFOLD;
  }else{
    return ((DBL_MAX_EXP + DBL_MANT_DIG - DIWIDTH + 1 + EXP_BIAS) - EXP(priX[0]))/DIWIDTH;
  }
  */
  return ((DBL_MAX_EXP + DBL_MANT_DIG - DIWIDTH + 1 + EXP_BIAS) - EXP(priX[0]))/DIWIDTH;
}
