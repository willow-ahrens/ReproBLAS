/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <math.h>
#include <float.h>

#include "../Common/Common.h"
#include "indexed.h"

/**
 * @internal
 * @brief  Add manually specified indexed double precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dmdmadd(const int fold, const double *manX, const int incmanX, const double *carX, const int inccarX, double* manY, const int incmanY, double* carY, const int inccarY) {
  int i;
  int shift;

  if (manX[0] == 0.0)
    return;

  if (manY[0] == 0.0) {
    for (i = 0; i < fold; i++) {
      manY[i*incmanY] = manX[i*incmanX];
      carY[i*inccarY] = carX[i*inccarX];
    }
    return;
  }

  if (isinf(manX[0]) || isnan(manX[0]) || isinf(manY[0]) || isnan(manY[0])) {
    manY[0] += manX[0];
    return;
  }

  shift = dmindex(manY) - dmindex(manX);
  if(shift > 0){
    //shift Y upwards and add X to Y
    for (i = fold - 1; i >= shift; i--) {
      manY[i*incmanY] = manX[i*incmanX] + (manY[(i - shift)*incmanY] - 1.5*ufp(manY[(i - shift)*incmanY]));
      carY[i*inccarY] = carX[i*inccarX] + carY[(i - shift)*inccarY];
    }
    for (i = 0; i < shift && i < fold; i++) {
      manY[i*incmanY] = manX[i*incmanX];
      carY[i*inccarY] = carX[i*inccarX];
    }
  }else{
    //shift X upwards and add X to Y
    for (i = 0 - shift; i < fold; i++) {
      manY[i*incmanY] += manX[(i + shift)*incmanX] - 1.5*ufp(manX[(i + shift)*incmanX]);
      carY[i*inccarY] += carX[(i + shift)*inccarX];
    }
  }

  dmrenorm(fold, manY, incmanY, carY, inccarY);
}

/**
 * @brief  Add indexed double precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void didiadd(const int fold, const double_indexed *X, double_indexed *Y){
  dmdmadd(fold, X, 1, X + fold, 1, Y, 1, Y + fold, 1);
}

/**
 * @internal
 * @brief  Add manually specified indexed complex double precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zmzmadd(const int fold, const double *manX, const int incmanX, const double *carX, const int inccarX, double* manY, const int incmanY, double* carY, const int inccarY) {
  dmdmadd(fold, manX, 2 * incmanX, carX, 2 * inccarX, manY, 2 * incmanY, carY, 2 * inccarY);
  dmdmadd(fold, manX + 1, 2 * incmanX, carX + 1, 2 * inccarX, manY + 1, 2 * incmanY, carY + 1, 2 * inccarY);
}

/**
 * @brief  Add indexed complex double precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void ziziadd(const int fold, const double_complex_indexed *X, double_complex_indexed *Y){
  zmzmadd(fold, X, 1, X + 2 * fold, 1, Y, 1, Y + 2 * fold, 1);
}

/**
 * @internal
 * @brief  Add double precision to suitably indexed manually specified indexed double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called dmdupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of dicapacity() elements into Y. After calling dmddeposit() on an indexed type, you must renormalize the indexed type with dmrenorm().
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dmddeposit(const int fold, const double X, double *manY, const int incmanY){
  double M;
  double x = X * dmcompression();
  long_double q;
  int i;

  if (isinf(x) || isnan(x) || isinf(manY[0]) || isnan(manY[0])) {
    manY[0] += x;
    return;
  }

  for (i = 0; i < fold - 1; i++) {
    M = manY[i * incmanY];
    q.d = x;
    q.l |= 1;
    q.d += M;
    manY[i * incmanY] = q.d;
    M -= q.d;
    x += M;
  }
  q.d = x;
  q.l |= 1;
  manY[i * incmanY] += q.d;
}

/**
 * @brief  Add double precision to suitably indexed indexed double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called didupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of dicapacity() elements into Y. After calling diddeposit() on an indexed type, you must renormalize the indexed type with direnorm().
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void diddeposit(const int fold, const double X, double_indexed *Y){
  dmddeposit(fold, X, Y, 1);
}

/**
 * @internal
 * @brief  Add complex double precision to suitably indexed manually specified indexed complex double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called zmzupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of dicapacity() elements into Y. After calling zmzdeposit() on an indexed type, you must renormalize the indexed type with zmrenorm().
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zmzdeposit(const int fold, const void *X, double *manY, const int incmanY){
  double MR, MI;
  long_double qR, qI;
  int i;
  double xR = ((double*)X)[0] * dmcompression();
  double xI = ((double*)X)[1] * dmcompression();

  if (isinf(xR) || isnan(xR) || isinf(manY[0]) || isnan(manY[0])) {
    manY[0] += xR;
    dmddeposit(fold, xI, manY + incmanY, 2 * incmanY);
    return;
  }
  if (isinf(xI) || isnan(xI) || isinf(manY[1]) || isnan(manY[1])) {
    manY[1] += xI;
    dmddeposit(fold, xR, manY, 2 * incmanY);
    return;
  }

  for (i = 0; i < fold - 1; i++) {
    MR = manY[i * 2 * incmanY];
    MI = manY[i * 2 * incmanY + 1];
    qR.d = xR;
    qI.d = xI;
    qR.l |= 1;
    qI.l |= 1;
    qR.d += MR;
    qI.d += MI;
    manY[i * 2 * incmanY] = qR.d;
    manY[i * 2 * incmanY + 1] = qI.d;
    MR -= qR.d;
    MI -= qI.d;
    xR += MR;
    xI += MI;
  }
  qR.d = xR;
  qI.d = xI;
  qR.l |= 1;
  qI.l |= 1;
  manY[i * 2 * incmanY] += qR.d;
  manY[i * 2 * incmanY + 1] += qI.d;
}

/**
 * @brief  Add complex double precision to suitably indexed manually specified indexed complex double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called zizupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of dicapacity() elements into Y. After calling zizdeposit() on an indexed type, you must renormalize the indexed type with zirenorm().
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zizdeposit(const int fold, const void *X, double_complex_indexed *Y){
  zmzdeposit(fold, X, Y, 1);
}

/**
 * @internal
 * @brief  Add double precision to manually specified indexed double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dmdadd(const int fold, const double X, double *manY, const int incmanY, double *carY, const int inccarY){
  dmdupdate(fold, fabs(X), manY, incmanY, carY, inccarY);
  dmddeposit(fold, X, manY, incmanY);
  dmrenorm(fold, manY, incmanY, carY, inccarY);
}

/**
 * @brief  Add double precision to indexed double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void didadd(const int fold, const double X, double_indexed *Y){
  dmdadd(fold, X, Y, 1, Y + fold, 1);
}

/**
 * @internal
 * @brief  Add complex double precision to manually specified indexed complex double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zmzadd(const int fold, const void *X, double *manY, const int incmanY, double *carY, const int inccarY){
  double aX[2];
  aX[0] = fabs(((double*)X)[0]);
  aX[1] = fabs(((double*)X)[1]);
  zmzupdate(fold, aX, manY, incmanY, carY, inccarY);
  zmzdeposit(fold, X, manY, incmanY);
  zmrenorm(fold, manY, incmanY, carY, inccarY);
}

/**
 * @brief  Add complex double precision to indexed complex double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zizadd(const int fold, const void *X, double_complex_indexed *Y){
  zmzadd(fold, X, Y, 1, Y + 2 * fold, 1);
}
