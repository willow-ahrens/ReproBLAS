#include <idxd.h>

/**
 * @brief  Add complex double precision to suitably indexed manually specified indexed complex double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called #zizupdate() on Y with the maximum absolute value of all future elements you wish to deposit in Y, you can call #zizdeposit() to deposit a maximum of #idxd_DIENDURANCE elements into Y before renormalizing Y with #zirenorm(). After any number of successive calls of #zizdeposit() on Y, you must renormalize Y with #zirenorm() before using any other function on Y.
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param priY Y's primary vector
 * @param incpriY stride within Y's primary vector (use every incpriY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   10 Jun 2015
 */
void zizdeposit(const int fold, const void *X, double_complex_indexed *Y){
  zmzdeposit(fold, X, Y, 1);
}
