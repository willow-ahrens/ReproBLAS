#include <math.h>

#include <idxd.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Get a reproducible double precision scale
 *
 * For any given X, return a reproducible scaling factor Y of the form
 *
 * 2^(#DIWIDTH * z) where z is an integer
 *
 * such that
 *
 * Y * 2^(-@c DBL_MANT_DIG - #DIWIDTH - 1) < X < Y * 2\^(#DIWIDTH + 2)
 *
 * @param X double precision number to be scaled
 * @return reproducible scaling factor
 *
 * @author Peter Ahrens
 * @date   19 Jun 2015
 */
double idxd_dscale(const double X){
    int e = EXP(X);
    e = e < DIWIDTH ? DIWIDTH : e;
    e -= (e - EXP_BIAS - 1) % DIWIDTH;
    e -= EXP_BIAS;
    return ldexp(0.5, e);
}
