/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "indexed.h"
#include "../Common/Common.h"

/**
 * @brief  Add manually specified indexed double precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param repX X's rep vector
 * @param increpX stride within X's rep vector (use every increpX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dmdmadd(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double* repY, const int increpY, double* carY, const int inccarY) {
  int i;
  int shift;

  if (repX[0] == 0.0)
    return;

  if (repY[0] == 0.0) {
    for (i = 0; i < fold; i++) {
      repY[i*increpY] = repX[i*increpX];
      carY[i*inccarY] = carX[i*inccarX];
    }
    return;
  }

  shift = dmindex(repY) - dmindex(repX);
  if(shift > 0){
    //shift Y upwards and add X to Y
    for (i = fold - 1; i >= shift; i--) {
      repY[i*increpY] = repX[i*increpX] + (repY[(i - shift)*increpY] - 1.5*ufp(repY[(i - shift)*increpY]));
      carY[i*inccarY] = carX[i*inccarX] + carY[(i - shift)*inccarY];
    }
    for (i = 0; i < shift && i < fold; i++) {
      repY[i*increpY] = repX[i*increpX];
      carY[i*inccarY] = carX[i*inccarX];
    }
  }else{
    //shift X upwards and add X to Y
    for (i = 0 - shift; i < fold; i++) {
      repY[i*increpY] += repX[(i + shift)*increpX] - 1.5*ufp(repX[(i + shift)*increpX]);
      carY[i*inccarY] += carX[(i + shift)*inccarX];
    }
  }

  dmrenorm(fold, repY, increpY, carY, inccarY);
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
 * @brief  Add manually specified indexed complex double precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param repX X's rep vector
 * @param increpX stride within X's rep vector (use every increpX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zmzmadd(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double* repY, const int increpY, double* carY, const int inccarY) {
  dmdmadd(fold, repX, 2 * increpX, carX, 2 * inccarX, repY, 2 * increpY, carY, 2 * inccarY);
  dmdmadd(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX, repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
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
 * @brief  Add double precision to suitably indexed manually specified indexed double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called dmdupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of dicapacity() elements into Y. After calling dmddeposit() on an indexed type, you must renormalize the indexed type with dmrenorm().
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dmddeposit(const int fold, const double X, double *repY, const int increpY){
  double M;
  double x = X;
  long_double q;
  int i;
  for (i = 0; i < fold - 1; i++) {
    M = repY[i * increpY];
    q.d = x;
    q.l |= 1;
    q.d += M;
    repY[i * increpY] = q.d;
    M -= q.d;
    x += M;
  }
  q.d = x;
  q.l |= 1;
  repY[i * increpY] += q.d;
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
 * @brief  Add complex double precision to suitably indexed manually specified indexed complex double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called zmzupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of dicapacity() elements into Y. After calling zmzdeposit() on an indexed type, you must renormalize the indexed type with zmrenorm().
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zmzdeposit(const int fold, const void *X, double *repY, const int increpY){
  double MR, MI;
  long_double qR, qI;
  int i;
  double xR = ((double*)X)[0];
  double xI = ((double*)X)[1];

  for (i = 0; i < fold - 1; i++) {
    MR = repY[i * 2 * increpY];
    MI = repY[i * 2 * increpY + 1];
    qR.d = xR;
    qI.d = xI;
    qR.l |= 1;
    qI.l |= 1;
    qR.d += MR;
    qI.d += MI;
    repY[i * 2 * increpY] = qR.d;
    repY[i * 2 * increpY + 1] = qI.d;
    MR -= qR.d;
    MI -= qI.d;
    xR += MR;
    xI += MI;
  }
  MR = repY[i * 2 * increpY];
  MI = repY[i * 2 * increpY + 1];
  qR.d = xR;
  qI.d = xI;
  qR.l |= 1;
  qI.l |= 1;
  qR.d += MR;
  qI.d += MI;
  repY[i * 2 * increpY] = qR.d;
  repY[i * 2 * increpY + 1] = qI.d;
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
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zizdeposit(const int fold, const void *X, double_complex_indexed *Y){
  zmzdeposit(fold, X, Y, 1);
}

/**
 * @brief  Add double precision to manually specified indexed double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dmdadd(const int fold, const double X, double *repY, const int increpY, double *carY, const int inccarY){
  dmdupdate(fold, fabs(X), repY, increpY, carY, inccarY);
  dmddeposit(fold, X, repY, increpY);
  dmrenorm(fold, repY, increpY, carY, inccarY);
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
 * @brief  Add complex double precision to manually specified indexed complex double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zmzadd(const int fold, const void *X, double *repY, const int increpY, double *carY, const int inccarY){
  double aX[2];
  aX[0] = fabs(((double*)X)[0]);
  aX[1] = fabs(((double*)X)[1]);
  zmzupdate(fold, aX, repY, increpY, carY, inccarY);
  zmzdeposit(fold, X, repY, increpY);
  zmrenorm(fold, repY, increpY, carY, inccarY);
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
