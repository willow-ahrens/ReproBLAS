#include <indexed.h>

/**
 * @internal
 * @brief rescale manually specified indexed complex single precision sum of squares
 *
 * Rescale an indexed complex single precision sum of squares Y to Y' such that Y / (scaleY * scaleY) == Y' / (X * X) and #smindex(Y) == #sindex(1.0)
 *
 * Note that Y is assumed to have an index at least the index of 1.0, and that X >= scaleY
 *
 * @param fold the fold of the indexed types
 * @param X Y's new scaleY (X == #sscale(Y) for some @c float Y) (X >= scaleY)
 * @param scaleY Y's current scaleY (scaleY == #sscale(Y) for some @c float Y) (X >= scaleY)
 * @param manY Y's mantissa vector (#smindex(Y) >= #sindex(1.0))
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Peter Ahrens
 * @date   19 Jun 2015
 */
void cmsrescale(const int fold, const float X, const float scaleY, float *manY, const int incmanY, float *carY, const int inccarY){
  int i;
  float rescaleY;

  if(X == scaleY || X == 0.0 || scaleY == 0.0){
    return;
  }

  rescaleY = X/scaleY;
  rescaleY *= rescaleY;
  for(i = 0; i < fold; i++){
    manY[i * 2 * incmanY] /= rescaleY;
    manY[i * 2 * incmanY + 1] /= rescaleY;
    if(manY[i * incmanY] == 0.0){
      cmsupdate(fold - i, 0.0, manY + i * 2 * incmanY, incmanY, carY + i * 2 * inccarY, inccarY);
      return;
    }
  }
}
