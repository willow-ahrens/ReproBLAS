#include <indexed.h>

/**
 * @internal
 * @brief rescale manually specified indexed complex double precision sum of squares
 *
 * Rescale an indexed complex double precision sum of squares X to X' such that X / (scaleX * scaleX) == X' / (newscaleX * newscaleX) and #dmindex(X) == #dindex(1.0)
 *
 * Note that X is assumed to have an index smaller than the index of 1.0, and that newscaleX >= scaleX
 *
 * @param newscaleX X's new scaleX (newscaleX == #dscale(Y) for some @c double Y) (newscaleX >= scaleX)
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector (#dmindex(X) >= #dindex(1.0))
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param scaleX X's current scaleX (scaleX == #dscale(Y) for some @c double Y) (newscaleX >= scaleX)
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
void zmdrescale(const int fold, const double newscaleX, double *manX, const int incmanX, double *carX, const int inccarX, const double scaleX){
  int i;
  double rescaleX;

  if(newscaleX == scaleX){
    return;
  }

  rescaleX = newscaleX/scaleX;
  rescaleX *= rescaleX;
  for(i = 0; i < fold; i++){
    manX[i * incmanX] /= rescaleX;
    manX[i * incmanX + 1] /= rescaleX;
    if(manX[i * incmanX] == 0.0){
      zmdupdate(fold - i, 0.0, manX + 2 * i * incmanX, incmanX, carX + 2 * i * inccarX, inccarX);
      break;
    }
  }
}
