#include <indexed.h>
#include <float.h>

/**
 * @internal
 * @brief Get a reproducible scale Y such that #dindex(X/Y) == #dindex(1.0)
 *
 * @param X the input number to be scaled
 * @return reproducible scaling factor (if X == 0.0, returns smallest valid scale)
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
double discale(const double X){
  int exp;
  frexp(X, &exp);
  if(X == 0.0){
    exp = DBL_MIN_EXP;
  }
  return ldexp(0.5, (DBL_MAX_EXP - ((DBL_MAX_EXP - exp) / diwidth()) * diwidth()));
}
