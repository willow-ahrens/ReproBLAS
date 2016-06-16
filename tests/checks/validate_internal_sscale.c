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
  return "Validate sscale internally";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  (void)incY;
  int i;
  float scale;
  float scale_mantissa;
  int scale_exponent;

  util_random_seed();

  //allocate vector
  float *X = util_svec_alloc(N*(FLT_MAX_EXP - FLT_MIN_EXP), incX);

  //check
  for (i = 0; i < (FLT_MAX_EXP - FLT_MIN_EXP) * N; i++) {
    X[i * incX] = ldexpf(0.5 + 0.5 * util_drand(), (i/N) + FLT_MIN_EXP);
  }
  for (i = 0; i < N * (FLT_MAX_EXP - FLT_MIN_EXP); i++) {
    scale = idxd_sscale(X[i * incX]);
    scale_mantissa = frexpf(scale, &scale_exponent) * 2.0;
    scale_exponent -= 1;
    if(scale_mantissa != 1.0 || scale_exponent % SIWIDTH != 0){
      printf("idxd_sscale(%g) = %g = %g * 2^%d is not of form 2^(SIWIDTH * integer)\n", X[i * incX], scale, scale_mantissa, scale_exponent);
      return 1;
    }
    if(ldexpf(0.5, -FLT_MANT_DIG - SIWIDTH) * scale >= X[i * incX]){
      printf("idxd_sscale(%g) * 2^(−FLT_MANT_DIG − SIWIDTH − 1) >= %g\n", X[i * incX], X[i * incX]);
      printf("%g * 2^(−FLT_MANT_DIG − SIWIDTH − 1) >= %g\n", scale, X[i * incX]);
      return 1;
    }
    if(X[i * incX] >= ldexpf(0.5, SIWIDTH + 3) * scale){
      printf("%g >= idxd_sscale(%g) * 2^(SIWIDTH + 2)\n", X[i * incX], X[i * incX]);
      printf("%g >= %g * 2^(SIWIDTH + 2)\n", X[i * incX], scale);
      return 1;
    }
  }
  free(X);
  return 0;
}
