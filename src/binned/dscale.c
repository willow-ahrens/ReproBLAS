#include <math.h>

#include <binned.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Get a reproducible double precision scale
 *
 * For any given X, return a reproducible scaling factor Y of the form
 *
 * 2^(#DBWIDTH * z) where z is an integer
 *
 * such that
 *
 * Y * 2^(-@c DBL_MANT_DIG - #DBWIDTH - 1) < X < Y * 2\^(#DBWIDTH + 2)
 *
 * @param X double precision number to be scaled
 * @return reproducible scaling factor
 *
 * @author Willow Ahrens
 * @date   19 Jun 2015
 */
double binned_dscale(const double X){
    int e = EXP(X);
    e = e < DBWIDTH ? DBWIDTH : e;
    e -= (e - EXP_BIAS - 1) % DBWIDTH;
    e -= EXP_BIAS;
    return ldexp(0.5, e);
}
