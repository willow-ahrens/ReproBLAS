#include <idxd.h>

/**
 * @brief  Add complex double precision to suitably indexed indexed complex double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called #idxd_zizupdate() on Y with the maximum absolute value of all future elements you wish to deposit in Y, you can call #idxd_zizdeposit() to deposit a maximum of #idxd_DIENDURANCE elements into Y before renormalizing Y with #idxd_zirenorm(). After any number of successive calls of #idxd_zizdeposit() on Y, you must renormalize Y with #idxd_zirenorm() before using any other function on Y.
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   10 Jun 2015
 */
void idxd_zizdeposit(const int fold, const void *X, double_complex_indexed *Y){
  idxd_zmzdeposit(fold, X, Y, 1);
}
