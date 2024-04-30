#include <math.h>

#include <binned.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Get a reproducible single precision scale
 *
 * For any given X, return a reproducible scaling factor Y of the form
 *
 * 2^(#SBWIDTH * z) where z is an integer
 *
 * such that
 *
 * Y * 2^(-@c FLT_MANT_DIG - #SBWIDTH - 1) < X < Y * 2\^(#SBWIDTH + 2)
 *
 * @param X single precision number to be scaled
 * @return reproducible scaling factor
 *
 * @author Willow Ahrens
 * @date   19 Jun 2015
 */
float binned_sscale(const float X){
    int e = EXPF(X);
    e = e < SBWIDTH ? SBWIDTH : e;
    e -= (e - EXPF_BIAS - 1) % SBWIDTH;
    e -= EXPF_BIAS;
    return ldexpf(0.5, e);
}
