#include <indexed.h>

/**
 * @internal
 * @brief Add manually specified indexed single precision sacled sums of squares (Y += X)
 *
 * Performs the operation Y += X, where X and Y represent scaled sums of squares.
 *
 * @param fold the fold of the indexed types
 * @param scaleX scale of X (scaleX == sscale(Z) for some @c float Z)
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param scaleY scale of Y (scaleY == sscale(Z) for some @c double Z)
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @return updated scale of Y
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
#include <stdio.h>
float smsmaddsq(const int fold, const float scaleX, const float *manX, const int incmanX, const float *carX, const int inccarX, const float scaleY, float* manY, const int incmanY, float* carY, const int inccarY) {
  if (scaleX > scaleY){
    smsrescale(fold, scaleX, scaleY, manY, incmanY, carY, inccarY);
    smsmadd(fold, manX, incmanX, carX, inccarX, manY, incmanY, carY, inccarY);
    return scaleX;
  }else if(scaleX == scaleY){
    smsmadd(fold, manX, incmanX, carX, inccarX, manY, incmanY, carY, inccarY);
    return scaleX;
  }else{
    float_indexed *tmp_X = sialloc(fold);
    smsmset(fold, manX, incmanX, carX, inccarX, tmp_X, 1, tmp_X + fold, 1);
    smsrescale(fold, scaleY, scaleX, tmp_X, 1, tmp_X + fold, 1);
    smsmadd(fold, tmp_X, 1, tmp_X + fold, 1, manY, incmanY, carY, inccarY);
    return scaleY;
  }
}
