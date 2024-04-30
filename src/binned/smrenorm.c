#include <math.h>

#include <binned.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Renormalize manually specified binned single precision
 *
 * Renormalization keeps the primary vector within the necessary bins by shifting over to the carry vector
 *
 * @param fold the fold of the binned types
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   23 Sep 2015
 */
void binned_smrenorm(const int fold, float* priX, const int incpriX, float* carX, const int inccarX) {
  /*
  //reference version
  int i;
  float M;
  float priX0 = priX[0];

  if(priX0 == 0.0 || ISNANINFF(priX0)){
    return;
  }

  for (i = 0; i < fold; i++, priX += incpriX, carX += inccarX) {
    priX0 = priX[0];

    M = UFPF(priX0);

    if (priX0 >= (M * 1.75)) {
      priX[0] -= M * 0.25;
      carX[0] += 1;
    }
    else if (priX0 < (M * 1.5)) {
      priX[0] += M * 0.25;
      carX[0] -= 1;
    }
  }
  */
  /*
  //vectorizeable version
  int i;
  int_float tmp_renorm, tmp_c;
  int tmp;

  for (i = 0; i < fold; i++, priX += incpriX, carX += inccarX) {
    tmp_renorm.f = priX[0];

    tmp_c.f = priX[0];
    tmp_c.i &= ((1ull << (FLT_MANT_DIG - 3)) | (1ul << (FLT_MANT_DIG - 2)));
    tmp_c.i <<= (33 - FLT_MANT_DIG);
    carX[0] -= 0.5 * tmp_c.f;

    tmp = tmp_renorm.i & (1ull << (FLT_MANT_DIG - 3));
    tmp <<= 1;
    tmp_renorm.i |= tmp;
    tmp_renorm.i &= ~(1ul << (FLT_MANT_DIG - 3));
    priX[0] = tmp_renorm.f;
  }
  */
  int i;
  int_float tmp_renorm;

  if(priX[0] == 0.0 || ISNANINFF(priX[0])){
    return;
  }

  for (i = 0; i < fold; i++, priX += incpriX, carX += inccarX) {
    tmp_renorm.f = priX[0];

    carX[0] += (int)((tmp_renorm.i >> (FLT_MANT_DIG - 3)) & 3) - 2;

    tmp_renorm.i &= ~(1ul << (FLT_MANT_DIG - 3));
    tmp_renorm.i |= 1ul << (FLT_MANT_DIG - 2);
    priX[0] = tmp_renorm.f;
  }
}
