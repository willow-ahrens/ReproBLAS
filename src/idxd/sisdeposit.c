#include <idxd.h>

/**
 * @brief  Add single precision to suitably indexed indexed single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called #idxd_sisupdate() on Y with the maximum absolute value of all future elements you wish to deposit in Y, you can call #idxd_sisdeposit() to deposit a maximum of #idxd_SIENDURANCE elements into Y before renormalizing Y with #idxd_sirenorm(). After any number of successive calls of #idxd_sisdeposit() on Y, you must renormalize Y with #idxd_sirenorm() before using any other function on Y.
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   10 Jun 2015
 */
void idxd_sisdeposit(const int fold, const float X, float_indexed *Y){
  idxd_smsdeposit(fold, X, Y, 1);
}
