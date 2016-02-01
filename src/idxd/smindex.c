#include <math.h>

#include <idxd.h>

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
 * @date   23 Sep 2015
 */
int idxd_smindex(const float *priX){
  /*
  //reference version
  int exp;

  if(priX[0] == 0.0){
    return (FLT_MAX_EXP - FLT_MIN_EXP)/SIWIDTH + idxd_SIMAXFOLD;
  }else{
    frexpf(priX[0], &exp);
    if(exp == FLT_MAX_EXP){
      return 0;
    }
    return (FLT_MAX_EXP + FLT_MANT_DIG - SIWIDTH + 1 - exp)/SIWIDTH;
  }
  */
  return ((FLT_MAX_EXP + FLT_MANT_DIG - SIWIDTH + 1 + EXPF_BIAS) - EXPF(priX[0]))/SIWIDTH;
}
