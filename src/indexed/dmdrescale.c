#include <indexed.h>

/**
 * @internal
 * @brief rescaleX manually specified indexed double precision sum of squares
 *
 * RescaleXs an indexed sum of squares X to X' such that X / (scaleX * scaleX) == X' / (sclX * sclX) and #dmindex(X) == #dindex(1.0)
 *
 * Note that X is assumed to have an index large enough to contain unity, and that sclX >= scaleX
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector (#dmindex(X) <= #dindex(1.0))
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param scaleX X's current scaleX (scaleX == #discaleX(Y) for some @c double Y)
 * @param sclX X's new scaleX (sclX == #discaleX(Y) for some @c double Y) (scaleX >= sclX)
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
void dmdrescaleX(const int fold, double *manX, const int incmanX, double *carX, const int inccarX, const double scaleX, const double sclX){
  int i;
  double rescaleX;

  if(sclX == scaleX){
    return;
  }

  rescaleX = sclX/scaleX;
  rescaleX *= rescaleX;
  for(i = 0; i < fold; i++){
    manX[i * incmanX] /= rescaleX;
    if(manX[i * incmanX] == 0.0){
      dmdupdate(fold - i, 0.0, manX + i * incmanX, incmanX, carX + i * inccarX, inccarX);
      break;
    }
  }
}
