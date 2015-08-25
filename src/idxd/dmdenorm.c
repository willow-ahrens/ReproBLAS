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
int idxd_dmdenorm(const int fold, const double *priX){
  return EXP(priX[0]) - DIWIDTH * fold < DBL_MIN_EXP + DBL_MANT_DIG + EXP_BIAS + DIWIDTH;
}
