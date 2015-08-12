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
 * @date   19 Jun 2015
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
  int exp = EXP(X);
  if(exp == 0){
    if(X == 0.0){
      return idxd_DIMAXINDEX;
    }else{
      frexp(X, &exp);
      return MIN((DBL_MAX_EXP - exp)/DIWIDTH, idxd_DIMAXINDEX);
    }
  }
  return ((DBL_MAX_EXP + EXP_BIAS) - exp)/DIWIDTH;
}
