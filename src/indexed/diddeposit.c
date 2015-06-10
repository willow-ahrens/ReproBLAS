#include <indexed.h>

/**
 * @brief  Add double precision to suitably indexed indexed double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called #didupdate() on Y with the maximum absolute value of all future elements you wish to deposit in Y, you can call #diddeposit() to deposit a maximum of #DIENDURANCE elements into Y before renormalizing Y with #direnorm(). After any number of successive calls of #diddeposit() on Y, you must renormalize Y with #direnorm() before using any other function on Y.
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   10 Jun 2015
 */
void diddeposit(const int fold, const double X, double_indexed *Y){
  dmddeposit(fold, X, Y, 1);
}
