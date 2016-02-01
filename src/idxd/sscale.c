#include <math.h>

#include <idxd.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Get a reproducible single precision scale
 *
 * For any given X, return a reproducible scaling factor Y of the form
 *
 * 2^(#SIWIDTH * z) where z is an integer
 *
 * such that
 *
 * Y * 2^(-@c FLT_MANT_DIG - #SIWIDTH - 1) < X < Y * 2\^(#SIWIDTH + 2)
 *
 * @param X single precision number to be scaled
 * @return reproducible scaling factor
 *
 * @author Peter Ahrens
 * @date   19 Jun 2015
 */
float idxd_sscale(const float X){
    int e = EXPF(X);
    e = e < SIWIDTH ? SIWIDTH : e;
    e -= (e - EXPF_BIAS - 1) % SIWIDTH;
    e -= EXPF_BIAS;
    return ldexpf(0.5, e);
}
