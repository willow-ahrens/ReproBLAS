#include <math.h>

#include <idxd.h>

#include "../common/common.h"

/**
 * @internal
 * @brief Renormalize manually specified indexed double precision
 *
 * Renormalization keeps the primary vector within the necessary bins by shifting over to the carry vector
 *
 * @param fold the fold of the indexed types
 * @param priX X's primary vector
 * @param incpriX stride within X's primary vector (use every incpriX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   23 Sep 2015
 */
void idxd_dmrenorm(const int fold, double* priX, const int incpriX, double* carX, const int inccarX) {
  /*
  //reference version
  int i;
  double M;
  double priX0 = priX[0];

  if(priX0 == 0.0 || ISNANINF(priX0)){
    return;
  }

  for (i = 0; i < fold; i++, priX += incpriX, carX += inccarX) {
    priX0 = priX[0];

    M = UFP(priX0);

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
  long_double tmp_renorm, tmp_c;
  long tmp;

  for (i = 0; i < fold; i++, priX += incpriX, carX += inccarX) {
    tmp_renorm.d = priX[0];
    tmp_c.d = priX[0];

    tmp_c.l &= ((1ull << (DBL_MANT_DIG - 3)) | (1ull << (DBL_MANT_DIG - 2)));
    tmp_c.l <<= (65 - DBL_MANT_DIG);
    carX[0] -= 0.5 * tmp_c.d;

    tmp = tmp_renorm.l & (1ull << (DBL_MANT_DIG - 3));
    tmp <<= 1;
    tmp_renorm.l |= tmp;
    tmp_renorm.l &= ~(1ull << (DBL_MANT_DIG - 3));
    priX[0] = tmp_renorm.d;
  }
  */
  int i;
  long_double tmp_renorm;

  if(priX[0] == 0.0 || ISNANINF(priX[0])){
    return;
  }

  for (i = 0; i < fold; i++, priX += incpriX, carX += inccarX) {
    tmp_renorm.d = priX[0];

    carX[0] += (int)((tmp_renorm.l >> (DBL_MANT_DIG - 3)) & 3) - 2;

    tmp_renorm.l &= ~(1ull << (DBL_MANT_DIG - 3));
    tmp_renorm.l |= 1ull << (DBL_MANT_DIG - 2);
    priX[0] = tmp_renorm.d;
  }
}
