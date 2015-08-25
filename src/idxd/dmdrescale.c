#include <idxd.h>

/**
 * @internal
 * @brief rescale manually specified indexed double precision sum of squares
 *
 * Rescale an indexed double precision sum of squares Y to Y' such that Y / (scaleY * scaleY) == Y' / (X * X)
 *
 * Note that Y is assumed to have an index at least the index of 1.0, and that X >= scaleY
 *
 * @param fold the fold of the indexed types
 * @param X Y's new scaleY (X == #idxd_dscale(Y) for some @c double Y) (X >= scaleY)
 * @param scaleY Y's current scaleY (scaleY == #idxd_dscale(Y) for some @c double Y) (X >= scaleY)
 * @param priY Y's primary vector (#idxd_dmindex(Y) >= #idxd_dindex(1.0))
 * @param incpriY stride within Y's primary vector (use every incpriY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Peter Ahrens
 * @date   19 Jun 2015
 */
void idxd_dmdrescale(const int fold, const double X, const double scaleY, double *priY, const int incpriY, double *carY, const int inccarY){
  int i;
  double rescaleY;

  if(X == scaleY || X == 0.0 || scaleY == 0.0){
    return;
  }

  rescaleY = X/scaleY;
  rescaleY *= rescaleY;
  for(i = 0; i < fold; i++){
    priY[i * incpriY] /= rescaleY;
    if(priY[i * incpriY] == 0.0){
      idxd_dmdupdate(fold - i, 0.0, priY + i * incpriY, incpriY, carY + i * inccarY, inccarY);
      return;
    }
  }
}
