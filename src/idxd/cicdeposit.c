#include <idxd.h>

/**
 * @brief  Add complex single precision to suitably indexed indexed complex single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called #idxd_cicupdate() on Y with the maximum absolute value of all future elements you wish to deposit in Y, you can call #idxd_cicdeposit() to deposit a maximum of #idxd_SIENDURANCE elements into Y before renormalizing Y with #idxd_cirenorm(). After any number of successive calls of #idxd_cicdeposit() on Y, you must renormalize Y with #idxd_cirenorm() before using any other function on Y.
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   10 Jun 2015
 */
void idxd_cicdeposit(const int fold, const void *X, float_complex_indexed *Y){
  idxd_cmcdeposit(fold, X, Y, 1);
}
