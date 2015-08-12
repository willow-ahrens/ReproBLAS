#include <indexed.h>

/**
 * @brief  Add single precision to suitably indexed indexed single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called #sisupdate() on Y with the maximum absolute value of all future elements you wish to deposit in Y, you can call #sisdeposit() to deposit a maximum of #idxd_SIENDURANCE elements into Y before renormalizing Y with #sirenorm(). After any number of successive calls of #sisdeposit() on Y, you must renormalize Y with #sirenorm() before using any other function on Y.
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   10 Jun 2015
 */
void sisdeposit(const int fold, const float X, float_indexed *Y){
  smsdeposit(fold, X, Y, 1);
}
