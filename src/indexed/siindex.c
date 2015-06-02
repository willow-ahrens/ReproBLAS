/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <math.h>
#include <float.h>

#include "indexed.h"
#include "../../config.h"
#include "../Common/Common.h"

#define FLT_BIN_DIG        13
static float bins[(FLT_MAX_EXP - FLT_MIN_EXP)/FLT_BIN_DIG + MAX_FOLD]; //initialized in bins_initialize
static int   bins_initialized  = 0;                                    //initialized in bins_initialize

/**
 * @brief Get indexed single precision bin width
 *
 * @return bin width (in bits)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int siwidth() {
  return FLT_BIN_DIG;
}

/**
 * @brief Get indexed single precision deposit capacity
 *
 * The number of deposits that can be performed before a renorm is necessary. This function applies also to indexed complex single precision.
 *
 * @return deposit capacity
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
int sicapacity() {
  return 1 << (FLT_MANT_DIG - FLT_BIN_DIG - 2);
}

/**
 * @internal
 * @brief Get indexed single precision compression factor
 *
 * This factor is used to scale down inputs before deposition
 *
 * @return compression factor
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
float smcompression() {
  return 1.0/(1 << (FLT_MANT_DIG - FLT_BIN_DIG + 1));
}

/**
 * @internal
 * @brief Get indexed single precision expansion factor
 *
 * This factor is used to scale up inputs after deposition
 *
 * @return expansion factor
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
float smexpansion() {
  return 1.0*(1 << (FLT_MANT_DIG - FLT_BIN_DIG + 1));
}

/**
 * @brief Get indexed single precision summation error bound
 *
 * This is a bound on the absolute error of a summation using indexed types
 *
 * @param fold the fold of the indexed types
 * @param N the number of single precision floating point summands
 * @param X the maximum absolute value of the summands
 * @return error bound
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   21 May 2015
 */
float sibound(const int fold, const int N, const float X) {
  return X * ldexpf(0.5, (1 - fold)*(FLT_BIN_DIG - 1) + 1) * N;
}

/**
 * @internal
 * @brief Get a reproducible single precision scale
 *
 * The scaling factor Y returned for given X is the smallest value that will fit in X's bin (The smallest representable value with the same index as X)
 *
 * Perhaps the most useful property of this number is that 1.0 <= X/Y < 2^#siwidth()
 *
 * @param X single precision number to be scaled
 * @return reproducible scaling factor (if X == 0.0, returns smallest valid scale)
 *
 * @author Peter Ahrens
 * @date   1 Jun 2015
 */
float sscale(const float X){
  return ldexpf(0.5, MAX((FLT_MAX_EXP - FLT_BIN_DIG + 1 - sindex(X) * FLT_BIN_DIG), FLT_MIN_EXP));
}

/**
 * @internal
 * @brief Get index of manually specified indexed single precision
 *
 * The index of an indexed type is the bin that it corresponds to. Higher indicies correspond to smaller bins.
 *
 * @param manX X's mantissa vector
 * @return X's index
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   19 May 2015
 */
int smindex(const float *manX){
  int exp;

  if(manX[0] == 0.0){
    return (FLT_MAX_EXP - FLT_MIN_EXP)/FLT_BIN_DIG + MAX_FOLD;
  }else{
    frexpf(manX[0], &exp);
    if(exp == FLT_MAX_EXP){
      return 0;
    }
    return (FLT_MAX_EXP + FLT_MANT_DIG - FLT_BIN_DIG + 1 - exp)/FLT_BIN_DIG;
  }
}

/**
 * @internal
 * @brief Check if index of manually specified indexed single precision is 0
 *
 * A quick check to determine if the index is 0
 *
 * @param manX X's mantissa vector
 * @return 1 if x has index 0, 0 otherwise.
 *
 * @author Peter Ahrens
 * @date   19 May 2015
 */
int smindex0(const float *manX){
  int exp;

  frexpf(manX[0], &exp);
  if(exp == FLT_MAX_EXP){
    return 1;
  }
  return 0;
}

/**
 * @brief Get index of single precision
 *
 * The index of a non-indexed type is the smallest index an indexed type would need to have to sum it reproducibly. Higher indicies correspond to smaller bins.
 *
 * @param X scalar X
 * @return X's index
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   19 May 2015
 */
int sindex(const float X){
  int exp;

  if(X == 0.0){
    return (FLT_MAX_EXP - FLT_MIN_EXP)/FLT_BIN_DIG;
  }else{
    frexpf(X, &exp);
    return (FLT_MAX_EXP - exp)/FLT_BIN_DIG;
  }
}

static void bins_initialize() {
  int index;

  if (bins_initialized) {
    return;
  }

  bins[0] = ldexpf(0.75, FLT_MAX_EXP);
  for(index = 1; index <= (FLT_MAX_EXP - FLT_MIN_EXP)/FLT_BIN_DIG; index++){
    bins[index] = ldexpf(0.75, (FLT_MAX_EXP + FLT_MANT_DIG - FLT_BIN_DIG + 1 - index * FLT_BIN_DIG));
  }
  for(; index < (FLT_MAX_EXP - FLT_MIN_EXP)/FLT_BIN_DIG + MAX_FOLD; index++){
    bins[index] = bins[index - 1];
  }

  bins_initialized = 1;
}

/**
 * @brief Get single precision bin corresponding to index
 *
 * @param X index
 * @return bin (bin)
 *
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
float sbin(const int X){
  bins_initialize();

  return bins[X];
}

/**
 * @internal
 * @brief Set manually specified indexed single precision bins
 *
 * Set the manually specified indexed single precision X to be empty with the given index
 *
 * @param X index
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void smbin(const int fold, const int X, float *manY, const int incmanY, float *carY, const int inccarY) {
  int i;

  bins_initialize();

  for (i = 0; i < fold; i++) {
    manY[i * incmanY] = bins[X + i];
    carY[i * inccarY] = 0.0;
  }
}
