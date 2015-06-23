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
 * @param priY Y's primary vector (#smindex(Y) >= #sindex(1.0))
 * @param incpriY stride within Y's primary vector (use every incpriY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Peter Ahrens
 * @date   19 Jun 2015
 */
void cmsrescale(const int fold, const float X, const float scaleY, float *priY, const int incpriY, float *carY, const int inccarY){
  int i;
  float rescaleY;

  if(X == scaleY || X == 0.0 || scaleY == 0.0){
    return;
  }

  rescaleY = X/scaleY;
  rescaleY *= rescaleY;
  for(i = 0; i < fold; i++){
    priY[i * 2 * incpriY] /= rescaleY;
    priY[i * 2 * incpriY + 1] /= rescaleY;
    if(priY[i * incpriY] == 0.0){
      cmsupdate(fold - i, 0.0, priY + i * 2 * incpriY, incpriY, carY + i * 2 * inccarY, inccarY);
      return;
    }
  }
}
