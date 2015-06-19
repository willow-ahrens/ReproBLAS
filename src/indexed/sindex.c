#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @brief Get index of single precision
 *
 * The index of a non-indexed type is the smallest index an indexed type would need to have to sum it reproducibly. Higher indicies correspond to smaller bins.
 *
 * @param X scalar X
 * @return X's index
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   19 Jun 2015
 */
int sindex(const float X){
  /*
  int exp;

  if(X == 0.0){
    return (FLT_MAX_EXP - FLT_MIN_EXP)/SIWIDTH;
  }else{
    frexpf(X, &exp);
    return (FLT_MAX_EXP - exp)/SIWIDTH;
  }
  */
  return ((FLT_MAX_EXP + EXPF_BIAS) - EXPF(X))/SIWIDTH;
}
