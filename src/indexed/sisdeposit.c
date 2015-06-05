#include <indexed.h>

/**
 * @brief  Add single precision to suitably indexed indexed single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called sisupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of siendurance() elements into Y. After calling sisdeposit() on an indexed type, you must renormalize the indexed type with sirenorm().
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void sisdeposit(const int fold, const float X, float_indexed *Y){
  smsdeposit(fold, X, Y, 1);
}
