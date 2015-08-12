#include <math.h>

#include <idxd.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Check if indexed type has denormal bits
 *
 * A quick check to determine if calculations involving X cannot be performed with "denormals are zero"
 *
 * @param priX X's primary vector
 * @return >0 if x has denormal bits, 0 otherwise.
 *
 * @author Peter Ahrens
 * @date   23 Jun 2015
 */
int smdenorm(const int fold, const float *priX){
  return EXPF(priX[0]) - fold * SIWIDTH < FLT_MIN_EXP + FLT_MANT_DIG + EXPF_BIAS + SIWIDTH;
}
