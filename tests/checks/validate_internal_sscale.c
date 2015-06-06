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
  return "Validate sscale internally";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  (void)incY;
  int i;
  float scale;
  float ratio;
  float bound;

  util_random_seed();

  //allocate vector
  float *X = util_svec_alloc(N*(FLT_MAX_EXP - FLT_MIN_EXP), incX);

  //check
  for (i = 0; i < (FLT_MAX_EXP - FLT_MIN_EXP) * N; i++) {
    X[i * incX] = ldexpf(0.5 + 0.5 * util_drand48(), (i/N) + FLT_MIN_EXP);
  }
  for (i = 0; i < N * (FLT_MAX_EXP - FLT_MIN_EXP); i++) {
    scale = sscale(X[i * incX]);
    bound = ldexpf(0.5, SIWIDTH + 1);
    ratio = X[i * incX] * (1.0/scale);
    if(ratio < 1.0){
      printf("%g * (1.0/sscale(%g) !>= 1.0\n", X[i * incX], X[i * incX]);
      printf("%g * (1.0/%g) !>= 1.0\n", X[i * incX], scale);
      printf("%g * (%g) !>= 1.0\n", X[i * incX], 1.0/scale);
      printf("%g !>= 1.0\n", ratio);
      return 1;
    }
    if(ratio >= bound){
      printf("%g * (1.0/sscale(%g)) !< 2^SIWIDTH\n", X[i * incX], X[i * incX]);
      printf("%g * (1.0/%g) !< %g\n", X[i * incX], scale, bound);
      printf("%g * (%g) !< %g\n", X[i * incX], 1.0/scale, bound);
      printf("%g !< %g\n", ratio, bound);
      return 1;
    }
  }
  if(sindex(sscale(0.0)) != sindex(0.0)){
    printf("sindex(sscale(0.0)) != sindex(0.0)\n");
    return 1;
  }
  return 0;
}
