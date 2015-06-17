#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @brief Get index of double precision
 *
 * The index of a non-indexed type is the smallest index an indexed type would need to have to sum it reproducibly. Higher indicies correspond to smaller bins.
 *
 * @param X scalar X
 * @return X's index
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   19 May 2015
 */
int dindex(const double X){
  /*
  int exp;

  if(X == 0.0){
    return (DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH;
  }else{
    frexp(X, &exp);
    return (DBL_MAX_EXP - exp)/DIWIDTH;
  }
  */
  return ((DBL_MAX_EXP + EXP_BIAS) - EXP(X))/DIWIDTH;
}
