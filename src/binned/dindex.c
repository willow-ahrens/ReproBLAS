#include <math.h>

#include <binned.h>

#include "../common/common.h"

/**
 * @brief Get index of double precision
 *
 * The index of a non-binned type is the smallest index an binned type would need to have to sum it reproducibly. Higher indicies correspond to smaller bins.
 *
 * @param X scalar X
 * @return X's index
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   19 Jun 2015
 */
int binned_dindex(const double X){
  /*
  //reference version
  int exp;

  if(X == 0.0){
    return (DBL_MAX_EXP - DBL_MIN_EXP)/DBWIDTH;
  }else{
    frexp(X, &exp);
    return (DBL_MAX_EXP - exp)/DBWIDTH;
  }
  */
  int exp = EXP(X);
  if(exp == 0){
    if(X == 0.0){
      return binned_DBMAXINDEX;
    }else{
      frexp(X, &exp);
      return MIN((DBL_MAX_EXP - exp)/DBWIDTH, binned_DBMAXINDEX);
    }
  }
  return ((DBL_MAX_EXP + EXP_BIAS) - exp)/DBWIDTH;
}
