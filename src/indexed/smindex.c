#include <math.h>

#include <indexed.h>

#include "../common/common.h"

#include "../../config.h"

/**
 * @internal
 * @brief Get index of manually specified indexed single precision
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
int smindex(const float *priX){
  /*
  int exp;

  if(priX[0] == 0.0){
    return (FLT_MAX_EXP - FLT_MIN_EXP)/SIWIDTH + SIMAXFOLD;
  }else{
    frexpf(priX[0], &exp);
    if(exp == FLT_MAX_EXP){
      return 0;
    }
    return (FLT_MAX_EXP + FLT_MANT_DIG - SIWIDTH + 1 - exp)/SIWIDTH;
  }
  */
  /*
  if(priX[0] == 0.0){
    return (FLT_MAX_EXP - FLT_MIN_EXP)/SIWIDTH + SIMAXFOLD;
  }else{
    return ((FLT_MAX_EXP + FLT_MANT_DIG - SIWIDTH + 1 + EXPF_BIAS) - EXPF(priX[0]))/SIWIDTH;
  }
  */
  return ((FLT_MAX_EXP + FLT_MANT_DIG - SIWIDTH + 1 + EXPF_BIAS) - EXPF(priX[0]))/SIWIDTH;
}
