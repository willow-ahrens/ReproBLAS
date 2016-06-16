#include <idxdBLAS.h>
#include <idxd.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "../common/test_vecvec_header.h"
#include "../common/test_util.h"


int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  (void)argc;
  (void)argv;
  return "Validate dscale internally";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  (void)incY;
  int i;
  double scale;
  double scale_mantissa;
  int scale_exponent;

  util_random_seed();

  //allocate vector
  double *X = util_dvec_alloc(N*(DBL_MAX_EXP - DBL_MIN_EXP), incX);

  //check
  for (i = 0; i < (DBL_MAX_EXP - DBL_MIN_EXP) * N; i++) {
    X[i * incX] = ldexp(0.5 + 0.5 * util_drand(), (i/N) + DBL_MIN_EXP);
  }
  for (i = 0; i < N * (DBL_MAX_EXP - DBL_MIN_EXP); i++) {
    scale = idxd_dscale(X[i * incX]);
    scale_mantissa = frexp(scale, &scale_exponent) * 2.0;
    scale_exponent -= 1;
    if(scale_mantissa != 1.0 || scale_exponent % DIWIDTH != 0){
      printf("idxd_dscale(%g) = %g = %g * 2^%d is not of form 2^(DIWIDTH * integer)\n", X[i * incX], scale, scale_mantissa, scale_exponent);
      return 1;
    }
    if(ldexp(0.5, -DBL_MANT_DIG - DIWIDTH) * scale >= X[i * incX]){
      printf("idxd_dscale(%g) * 2^(−DBL_MANT_DIG − DIWIDTH − 1) >= %g\n", X[i * incX], X[i * incX]);
      printf("%g * 2^(−DBL_MANT_DIG − DIWIDTH − 1) >= %g\n", scale, X[i * incX]);
      return 1;
    }
    if(X[i * incX] >= ldexp(0.5, DIWIDTH + 3) * scale){
      printf("%g >= idxd_dscale(%g) * 2^(DIWIDTH + 2)\n", X[i * incX], X[i * incX]);
      printf("%g >= %g * 2^(DIWIDTH + 2)\n", X[i * incX], scale);
      return 1;
    }
  }
  free(X);
  return 0;
}
