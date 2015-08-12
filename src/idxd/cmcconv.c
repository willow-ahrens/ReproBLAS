#include <idxd.h>

/**
 * @internal
 * @brief Convert complex single precision to manually specified indexed complex single precision (X -> Y)
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param priY Y's primary vector
 * @param incpriY stride within Y's primary vector (use every incpriY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void idxd_cmcconv(const int fold, const void *X, float *priY, const int incpriY, float *carY, const int inccarY) {
  idxd_smsconv(fold, ((float*)X)[0], priY, incpriY * 2, carY, inccarY * 2);
  idxd_smsconv(fold, ((float*)X)[1], priY + 1, incpriY * 2, carY + 1, inccarY * 2);
}
