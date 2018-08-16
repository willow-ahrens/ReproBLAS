#include <binned.h>

/**
 * @brief  Add complex double precision to suitably binned binned complex double precision (Y += X)
 *
 * Performs the operation Y += X on an binned type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called #binned_zbzupdate() on Y with the maximum absolute value of all future elements you wish to deposit in Y, you can call #binned_zbzdeposit() to deposit a maximum of #binned_DBENDURANCE elements into Y before renormalizing Y with #binned_zbrenorm(). After any number of successive calls of #binned_zbzdeposit() on Y, you must renormalize Y with #binned_zbrenorm() before using any other function on Y.
 *
 * @param fold the fold of the binned types
 * @param X scalar X
 * @param Y binned scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   10 Jun 2015
 */
void binned_zbzdeposit(const int fold, const void *X, double_complex_binned *Y){
  binned_zmzdeposit(fold, X, Y, 1);
}
