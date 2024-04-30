#include <binned.h>

/**
 * @brief  Add single precision to suitably binned binned single precision (Y += X)
 *
 * Performs the operation Y += X on an binned type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called #binned_sbsupdate() on Y with the maximum absolute value of all future elements you wish to deposit in Y, you can call #binned_sbsdeposit() to deposit a maximum of #binned_SBENDURANCE elements into Y before renormalizing Y with #binned_sbrenorm(). After any number of successive calls of #binned_sbsdeposit() on Y, you must renormalize Y with #binned_sbrenorm() before using any other function on Y.
 *
 * @param fold the fold of the binned types
 * @param X scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   10 Jun 2015
 */
void binned_sbsdeposit(const int fold, const float X, float_binned *Y){
  binned_smsdeposit(fold, X, Y, 1);
}
