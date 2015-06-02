#include <indexed.h>

/**
 * @internal
 * @brief rescale manually specified indexed double precision sum of squares
 *
 * Rescale an indexed double precision sum of squares Y to Y' such that Y / (scaleY * scaleY) == Y' / (X * X) and #dmindex(Y) == #dindex(1.0)
 *
 * Note that Y is assumed to have an index smaller than the index of 1.0, and that X >= scaleY
 *
 * @param fold the fold of the indexed types
 * @param X Y's new scaleY (X == #dscale(Y) for some @c double Y) (X >= scaleY)
 * @param scaleY Y's current scaleY (scaleY == #dscale(Y) for some @c double Y) (X >= scaleY)
 * @param manY Y's mantissa vector (#dmindex(Y) >= #dindex(1.0))
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
void dmdrescale(const int fold, const double X, const double scaleY, double *manY, const int incmanY, double *carY, const int inccarY){
  int i;
  double rescaleY;

  if(X == scaleY || X == 0.0 || scaleY == 0.0){
    return;
  }

  rescaleY = X/scaleY;
  rescaleY *= rescaleY;
  for(i = 0; i < fold; i++){
    manY[i * incmanY] /= rescaleY;
    if(manY[i * incmanY] == 0.0){
      dmdupdate(fold - i, 0.0, manY + i * incmanY, incmanY, carY + i * inccarY, inccarY);
      break;
    }
  }
}
