#include <indexed.h>

/**
 * @internal
 * @brief rescale manually specified indexed complex single precision sum of squares
 *
 * Rescale an indexed complex single precision sum of squares X to X' such that X / (scaleX * scaleX) == X' / (newscaleX * newscaleX) and #smindex(X) == #sindex(1.0)
 *
 * Note that X is assumed to have an index smaller than the index of 1.0, and that newscaleX >= scaleX
 *
 * @param newscaleX X's new scaleX (newscaleX == #sscale(Y) for some @c float Y) (newscaleX >= scaleX)
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector (#smindex(X) >= #sindex(1.0))
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param scaleX X's current scaleX (scaleX == #sscale(Y) for some @c float Y) (newscaleX >= scaleX)
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
void cmsrescale(const int fold, const float newscaleX, float *manX, const int incmanX, float *carX, const int inccarX, const float scaleX){
  int i;
  float rescaleX;

  if(newscaleX == scaleX){
    return;
  }

  rescaleX = newscaleX/scaleX;
  rescaleX *= rescaleX;
  for(i = 0; i < fold; i++){
    manX[i * incmanX] /= rescaleX;
    manX[i * incmanX + 1] /= rescaleX;
    if(manX[i * incmanX] == 0.0){
      cmsupdate(fold - i, 0.0, manX + 2 * i * incmanX, incmanX, carX + 2 * i * inccarX, inccarX);
      break;
    }
  }
}
