#include <binned.h>

/**
 * @internal
 * @brief Add manually specified binned single precision scaled sums of squares (Y += X)
 *
 * Performs the operation Y += X, where X and Y represent scaled sums of squares.
 *
 * @param fold the fold of the binned types
 * @param scaleX scale of X (scaleX == #binned_sscale (Z) for some @c float Z)
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param scaleY scale of Y (scaleY == #binned_sscale (Z) for some @c double Z)
 * @param priY Y's primary vector
 * @param incpriY stride within Y's primary vector (use every incpriY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @return updated scale of Y
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
#include <stdio.h>
float binned_smsmaddsq(const int fold, const float scaleX, const float *priX, const int incpriX, const float *carX, const int inccarX, const float scaleY, float *priY, const int incpriY, float *carY, const int inccarY) {
  if (scaleX > scaleY){
    binned_smsrescale(fold, scaleX, scaleY, priY, incpriY, carY, inccarY);
    binned_smsmadd(fold, priX, incpriX, carX, inccarX, priY, incpriY, carY, inccarY);
    return scaleX;
  }else if(scaleX == scaleY){
    binned_smsmadd(fold, priX, incpriX, carX, inccarX, priY, incpriY, carY, inccarY);
    return scaleX;
  }else{
    float_binned *tmp_X = binned_sballoc(fold);
    binned_smsmset(fold, priX, incpriX, carX, inccarX, tmp_X, 1, tmp_X + fold, 1);
    binned_smsrescale(fold, scaleY, scaleX, tmp_X, 1, tmp_X + fold, 1);
    binned_smsmadd(fold, tmp_X, 1, tmp_X + fold, 1, priY, incpriY, carY, inccarY);
    return scaleY;
  }
}
