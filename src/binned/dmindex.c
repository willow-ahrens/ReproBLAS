#include <math.h>

#include <binned.h>

#include "../common/common.h"

#include "../../config.h"

/**
 * @internal
 * @brief Get index of manually specified binned double precision
 *
 * The index of an binned type is the bin that it corresponds to. Higher indicies correspond to smaller bins.
 *
 * @param priX X's primary vector
 * @return X's index
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   23 Sep 2015
 */
int binned_dmindex(const double *priX){
  /*
  //reference version
  int exp;

  if(priX[0] == 0.0){
    return (DBL_MAX_EXP - DBL_MIN_EXP)/DBWIDTH + binned_DBMAXFOLD;
  }else{
    frexp(priX[0], &exp);
    if(exp == DBL_MAX_EXP){
      return 0;
    }
    return (DBL_MAX_EXP + DBL_MANT_DIG - DBWIDTH + 1 - exp)/DBWIDTH;
  }
  */
  return ((DBL_MAX_EXP + DBL_MANT_DIG - DBWIDTH + 1 + EXP_BIAS) - EXP(priX[0]))/DBWIDTH;
}
