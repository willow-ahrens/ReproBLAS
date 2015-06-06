#include <indexedBLAS.h>
#include <indexed.h>
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
  return "Verify sindex";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  (void)incY;
  int i;
  int index;

  util_random_seed();

  //allocate vector
  float *X = util_svec_alloc(N*(FLT_MAX_EXP - FLT_MIN_EXP) + 1, incX);

  //check
  for (i = 0; i < (FLT_MAX_EXP - FLT_MIN_EXP) * N; i++) {
    X[i * incX] = ldexpf(0.5 + 0.5 * util_drand48(), (i/N) + FLT_MIN_EXP);
  }
  X[i * incX] = 0.0;
  for (i = 0; i < N * (FLT_MAX_EXP - FLT_MIN_EXP) + 1; i++) {
    index = sindex(X[i * incX]);
    if (index == 0){
      X[i * incX] *= SMCOMPRESSION;
    }
    if (X[i * incX] != 0.0 && sbin(index) / ldexpf(0.75, FLT_MANT_DIG) > 2 * fabsf(X[i * incX])){
      printf("2 * |X| !>= 2^(i * W)\n");
      printf("2 * %g !>= %g\n", fabsf(X[i * incX]), sbin(index + 1)/ldexpf(0.75, FLT_MANT_DIG - SIWIDTH));
      return 1;
    }
    if (sbin(index) / ldexpf(0.75, FLT_MANT_DIG - SIWIDTH) <= 2 * fabsf(X[i * incX])){
      printf("2 * |X| !< 2^((i + 1) * W)\n");
      printf("2 * %g !< %g\n", fabsf(X[i * incX]) * SMCOMPRESSION, sbin(index)/ldexpf(0.75, FLT_MANT_DIG - SIWIDTH));
      return 1;
    }
  }
  return 0;
}
