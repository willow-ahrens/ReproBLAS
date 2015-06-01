#include <indexed.h>
#include <float.h>

/**
 * @internal
 * @brief Get a reproducible single precision scale Y representing the magnitude of elements in X's bin.
 *
 * A property of this value is that #sindex(X/Y) == #sindex(1.0)
 *
 * @param X single precision number to be scaled
 * @return reproducible scaling factor (if X == 0.0, returns smallest valid scale)
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
float sscale(const float X){
  int exp;
  frexpf(X, &exp);
  if(X == 0.0){
    exp = FLT_MIN_EXP;
  }
  return ldexpf(0.5, (FLT_MAX_EXP - ((FLT_MAX_EXP - exp) / siwidth()) * siwidth()));
}
