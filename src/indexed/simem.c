#include <stdlib.h>
#include <string.h>

#include "indexed.h"

/**
 * @brief indexed single precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in bytes) of the indexed type
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
size_t sisize(const int fold){
  return 2*fold*sizeof(float);
}

/**
 * @brief indexed complex single precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in bytes) of the indexed type
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
size_t cisize(const int fold){
  return 4*fold*sizeof(float);
}

/**
 * @brief indexed single precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in @c float) of the indexed type
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int sinum(const int fold){
  return 2*fold;
}

/**
 * @brief indexed complex single precision size
 *
 * @param fold the fold of the indexed type
 * @return the size (in @c float) of the indexed type
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int cinum(const int fold){
  return 4*fold;
}

/**
 * @brief indexed single precision allocation
 *
 * @param fold the fold of the indexed type
 * @return a freshly allocated indexed type. (free with @c free())
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
float_indexed *sialloc(const int fold){
  return (float_indexed*)malloc(sisize(fold));
}

/**
 * @brief indexed complex single precision allocation
 *
 * @param fold the fold of the indexed type
 * @return a freshly allocated indexed type. (free with @c free())
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
float_complex_indexed *cialloc(const int fold){
  return (float_complex_indexed*)malloc(cisize(fold));
}

/**
 * @internal
 * @brief Set manually specified indexed single precision (Y = X)
 *
 * Performs the operation Y = X
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
void smsmset(const int fold, const float *manX, const int incmanX, const float *carX, const int inccarX, float *manY, const int incmanY, float *carY, const int inccarY){
  int i;
  for(i = 0; i < fold; i++){
    manY[i * incmanY] = manX[i * incmanX];
    carY[i * inccarY] = carX[i * inccarX];
  }
}

/**
 * @brief Set indexed single precision (Y = X)
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
void sisiset(const int fold, const float_indexed *X, float_indexed *Y){
  memcpy(Y, X, sisize(fold));
}

/**
 * @internal
 * @brief Set manually specified indexed complex single precision (Y = X)
 *
 * Performs the operation Y = X
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
void cmcmset(const int fold, const float *manX, const int incmanX, const float *carX, const int inccarX, float *manY, const int incmanY, float *carY, const int inccarY){
  smsmset(fold, manX, 2 * incmanX, carX, 2 * inccarX, manY, 2 * incmanY, carY, 2 * inccarY);
  smsmset(fold, manX + 1, 2 * incmanX, carX + 1, 2 * inccarX, manY + 1, 2 * incmanY, carY + 1, 2 * inccarY);
}

/**
 * @brief Set indexed complex single precision (Y = X)
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
void ciciset(const int fold, const float_complex_indexed *X, float_complex_indexed *Y){
  memcpy(Y, X, cisize(fold));
}

/**
 * @internal
 * @brief Set manually specified indexed complex single precision to manually specified indexed single precision (Y = X)
 *
 * Performs the operation Y = X
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
void cmsmset(const int fold, const float *manX, const int incmanX, const float *carX, const int inccarX, float *manY, const int incmanY, float *carY, const int inccarY){
  smsmset(fold, manX, incmanX, carX, inccarX, manY, 2 * incmanY, carY, 2 * inccarY);
  smsetzero(fold, manY + 1, 2 * incmanY, carY + 1, 2 * inccarY);
}

/**
 * @brief Set indexed complex single precision to indexed single precision (Y = X)
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
void cisiset(const int fold, const float_indexed *X, float_complex_indexed *Y){
  cmsmset(fold, X, 1, X + fold, 1, Y, 1, Y + 2 * fold, 1);
}

/**
 * @internal
 * @brief Set manually specified indexed single precision to 0 (X = 0)
 *
 * Performs the operation X = 0
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void smsetzero(const int fold, float *manX, const int incmanX, float *carX, const int inccarX){
  int i;
  for(i = 0; i < fold; i++){
    manX[i * incmanX] = 0.0;
    carX[i * inccarX] = 0.0;
  }
}

/**
 * @brief Set indexed single precision to 0 (X = 0)
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
void sisetzero(const int fold, float_indexed *X){
  memset(X, 0, sisize(fold));
}

/**
 * @internal
 * @brief Set manually specified indexed complex single precision to 0 (X = 0)
 *
 * Performs the operation X = 0
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmsetzero(const int fold, float *manX, const int incmanX, float *carX, const int inccarX){
  smsetzero(fold, manX, 2 * incmanX, carX, 2 * inccarX);
  smsetzero(fold, manX + 1, 2 * incmanX, carX + 1, 2 * inccarX);
}

/**
 * @brief Set indexed single precision to 0 (X = 0)
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
void cisetzero(const int fold, float_complex_indexed *X){
  memset(X, 0, cisize(fold));
}
