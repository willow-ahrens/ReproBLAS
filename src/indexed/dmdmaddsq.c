#include <indexed.h>

/**
 * @internal
 * @brief Add manually specified indexed double precision sacled sums of squares (Y += X)
 *
 * Performs the operation Y += X, where X and Y represent scaled sums of squares.
 *
 * @param fold the fold of the indexed types
 * @param scaleX scale of X (scaleX == dscale(Z) for some @c double Z)
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param scaleY scale of Y (scaleY == dscale(Z) for some @c double Z)
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
double dmdmaddsq(const int fold, const double scaleX, const double *manX, const int incmanX, const double *carX, const int inccarX, const double scaleY, double* manY, const int incmanY, double* carY, const int inccarY) {
  if (scaleX > scaleY){
    dmdrescale(fold, scaleX, manY, incmanY, carY, inccarY, scaleY);
    dmdmadd(fold, manX, incmanX, carX, inccarX, manY, incmanY, carY, inccarY);
    return scaleX;
  }else if(scaleX == scaleY){
    dmdmadd(fold, manX, incmanX, carX, inccarX, manY, incmanY, carY, inccarY);
    return scaleX;
  }else{
    double_indexed tmp_X = dialloc(fold);
    dmdmset(fold, manX, incmanX, carX, inccarX, tmp_X, 1, tmp_X + fold, 1);
    dmdrescale(fold, scaleY, tmpX, 1, tmpX + fold, 1, scaleX);
    dmdmadd(fold, tmp_X, 1, tmp_X + fold, 1, manY, incmanY, carY, inccarY);
    return scaleY;
  }
}
