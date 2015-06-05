#include <indexed.h>

/**
 * @brief  Add complex double precision to suitably indexed manually specified indexed complex double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called zizupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of diendurance() elements into Y. After calling zizdeposit() on an indexed type, you must renormalize the indexed type with zirenorm().
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zizdeposit(const int fold, const void *X, double_complex_indexed *Y){
  zmzdeposit(fold, X, Y, 1);
}
