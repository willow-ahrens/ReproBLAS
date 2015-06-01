#include <indexed.h>

/**
 * @internal
 * @brief rescale manually specified indexed double precision
 *
 * Rescales X to X' such that X / (ScaleX * ScaleX) == X' / (NewScaleX * NewScaleX) and #dmindex(X) == #dindex(1.0)
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector (#dmindex(X) <= #dindex(1.0))
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param ScaleX X's current scale (ScaleX == #discale(Y) for some @c double Y)
 * @param NewScaleX X's new scale (NewScaleX == #discale(Y) for some @c double Y) (NewScaleX >= ScaleX)
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
void dmdrescale(const int fold, double *manX, const int incmanX, double *carX, const int inccarX, const double ScaleX, const double NewScaleX){
  int i;
  double ReScaleX;

  if(ScaleX == NewScaleX){
    return;
  }

  ReScaleX = NewScaleX/ScaleX;
  ReScaleX *= ReScaleX;
  for(i = 0; i < fold; i++){
    manX[i * incmanX] /= ReScaleX;
    if(manX[i * incmanX] == 0.0){
      dmdupdate(fold - i, 0.0, manX + i * incmanX, incmanX, carX + i * inccarX, inccarX);
      break;
    }
  }
}
