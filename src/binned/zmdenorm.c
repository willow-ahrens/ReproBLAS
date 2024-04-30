#include <math.h>

#include <binned.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Check if binned type has denormal bits
 *
 * A quick check to determine if calculations involving X cannot be performed with "denormals are zero"
 *
 * @param fold the fold of the binned type
 * @param priX X's primary vector
 * @return >0 if x has denormal bits, 0 otherwise.
 *
 * @author Willow Ahrens
 * @date   23 Jun 2015
 */
int binned_zmdenorm(const int fold, const double *priX){
  return binned_dmdenorm(fold, priX) || binned_dmdenorm(fold, priX + 1);
}
