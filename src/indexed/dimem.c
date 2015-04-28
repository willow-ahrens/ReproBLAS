#include <stdlib.h>
#include <string.h>

#include "indexed.h"

/**
 * @brief indexed double precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in bytes) of the indexed type
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
size_t disize(const int fold){
  return 2*fold*sizeof(double);
}

/**
 * @brief indexed complex double precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in bytes) of the indexed type
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
size_t zisize(const int fold){
  return 4*fold*sizeof(double);
}

/**
 * @brief indexed double precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in @c double) of the indexed type
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int dinum(const int fold){
  return 2*fold;
}

/**
 * @brief indexed complex double precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in @c double) of the indexed type
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int zinum(const int fold){
  return 4*fold;
}

/**
 * @brief indexed double precision allocation
 *
 * @param fold the fold of the indexed type
 * @return a freshly allocated indexed type. (free with @c free())
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
double_indexed *dialloc(const int fold){
  return (double_indexed*)malloc(disize(fold));
}

/**
 * @brief indexed complex double precision allocation
 *
 * @param fold the fold of the indexed type
 * @return a freshly allocated indexed type. (free with @c free())
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
double_complex_indexed *zialloc(const int fold){
  return (double_complex_indexed*)malloc(zisize(fold));
}

/**
 * @internal
 * @brief Set manually specified indexed double precision (Y = X)
 *
 * Performs the operation Y = X
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
void dmdmset(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double *repY, const int increpY, double *carY, const int inccarY){
  int i;
  for(i = 0; i < fold; i++){
    repY[i * increpY] = repX[i * increpX];
    carY[i * inccarY] = carX[i * inccarX];
  }
}

/**
 * @brief Set indexed double precision (Y = X)
 *
 * Performs the operation Y = X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void didiset(const int fold, const double_indexed *X, double_indexed *Y){
  memcpy(Y, X, disize(fold));
}

/**
 * @internal
 * @brief Set manually specified indexed complex double precision (Y = X)
 *
 * Performs the operation Y = X
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
void zmzmset(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double *repY, const int increpY, double *carY, const int inccarY){
  dmdmset(fold, repX, 2 * increpX, carX, 2 * inccarX, repY, 2 * increpY, carY, 2 * inccarY);
  dmdmset(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX, repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}

/**
 * @brief Set indexed complex double precision (Y = X)
 *
 * Performs the operation Y = X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void ziziset(const int fold, const double_complex_indexed *X, double_complex_indexed *Y){
  memcpy(Y, X, zisize(fold));
}

/**
 * @internal
 * @brief Set manually specified indexed complex double precision to manually specified indexed double precision (Y = X)
 *
 * Performs the operation Y = X
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
void zmdmset(const int fold, const double *repX, const int increpX, const double *carX, const int inccarX, double *repY, const int increpY, double *carY, const int inccarY){
  dmdmset(fold, repX, increpX, carX, inccarX, repY, 2 * increpY, carY, 2 * inccarY);
  dmsetzero(fold, repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}


/**
 * @brief Set indexed complex double precision to indexed double precision (Y = X)
 *
 * Performs the operation Y = X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zidiset(const int fold, const double_indexed *X, double_complex_indexed *Y){
  zmdmset(fold, X, 1, X + fold, 1, Y, 1, Y + 2 * fold, 1);
}

/**
 * @internal
 * @brief Set manually specified indexed double precision to 0 (X = 0)
 *
 * Performs the operation X = 0
 *
 * @param fold the fold of the indexed types
 * @param repX X's rep vector
 * @param increpX stride within X's rep vector (use every increpX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dmsetzero(const int fold, double *repX, const int increpX, double *carX, const int inccarX){
  int i;
  for(i = 0; i < fold; i++){
    repX[i * increpX] = 0.0;
    carX[i * inccarX] = 0.0;
  }
}

/**
 * @brief Set indexed double precision to 0 (X = 0)
 *
 * Performs the operation X = 0
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void disetzero(const int fold, double_indexed *X){
  memset(X, 0, disize(fold));
}

/**
 * @internal
 * @brief Set manually specified indexed complex double precision to 0 (X = 0)
 *
 * Performs the operation X = 0
 *
 * @param fold the fold of the indexed types
 * @param repX X's rep vector
 * @param increpX stride within X's rep vector (use every increpX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zmsetzero(const int fold, double *repX, const int increpX, double *carX, const int inccarX){
  dmsetzero(fold, repX, 2 * increpX, carX, 2 * inccarX);
  dmsetzero(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX);
}

/**
 * @brief Set indexed double precision to 0 (X = 0)
 *
 * Performs the operation X = 0
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void zisetzero(const int fold, double_complex_indexed *X){
  memset(X, 0, zisize(fold));
}
