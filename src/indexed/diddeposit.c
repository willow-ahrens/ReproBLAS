#include <indexed.h>

/**
 * @brief  Add double precision to suitably indexed indexed double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called didupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of DIENDURANCE elements into Y. After calling diddeposit() on an indexed type, you must renormalize the indexed type with direnorm().
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void diddeposit(const int fold, const double X, double_indexed *Y){
  dmddeposit(fold, X, Y, 1);
}
