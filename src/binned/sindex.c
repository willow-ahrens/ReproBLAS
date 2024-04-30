#include <math.h>

#include <binned.h>

#include "../common/common.h"

/**
 * @brief Get index of single precision
 *
 * The index of a non-binned type is the smallest index an binned type would need to have to sum it reproducibly. Higher indicies correspond to smaller bins.
 *
 * @param X scalar X
 * @return X's index
 *
 * @author Willow Ahrens
 * @author Hong Diep Nguyen
 * @date   19 Jun 2015
 */
int binned_sindex(const float X){
  /*
  //reference version
  int exp;

  if(X == 0.0){
    return (FLT_MAX_EXP - FLT_MIN_EXP)/SBWIDTH;
  }else{
    frexpf(X, &exp);
    return (FLT_MAX_EXP - exp)/SBWIDTH;
  }
  */
  int exp = EXPF(X);
  if(exp == 0){
    if(X == 0.0){
      return binned_SBMAXINDEX;
    }else{
      frexpf(X, &exp);
      return MIN((FLT_MAX_EXP - exp)/SBWIDTH, binned_SBMAXINDEX);
    }
  }
  return ((FLT_MAX_EXP + EXPF_BIAS) - exp)/SBWIDTH;
}
