#include <math.h>

#include <binned.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Check if index of manually specified binned double precision is 0
 *
 * A quick check to determine if the index is 0
 *
 * @param priX X's primary vector
 * @return >0 if x has index 0, 0 otherwise.
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
int binned_dmindex0(const double *priX){
  /*
  //reference version
  int exp;

  frexp(priX[0], &exp);
  if(exp == DBL_MAX_EXP){
    return 1;
  }
  return 0;
  */
  return EXP(priX[0]) == DBL_MAX_EXP + EXP_BIAS;
}
