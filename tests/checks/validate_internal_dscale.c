#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

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
  double bound;
  double ratio;

  util_random_seed();

  //allocate vector
  double *X = util_dvec_alloc(N*(DBL_MAX_EXP - DBL_MIN_EXP), incX);

  //check
  for (i = 0; i < (DBL_MAX_EXP - DBL_MIN_EXP) * N; i++) {
    X[i * incX] = ldexp(0.5 + 0.5 * util_drand48(), (i/N) + DBL_MIN_EXP);
  }
  for (i = 0; i < N * (DBL_MAX_EXP - DBL_MIN_EXP); i++) {
    scale = dscale(X[i * incX]);
    bound = ldexp(0.5, diwidth() + 1);
    ratio = X[i * incX] * (1.0/scale);
    if(ratio < 1.0){
      printf("%g * (1.0/dscale(%g) !>= 1.0\n", X[i * incX], X[i * incX]);
      printf("%g * (1.0/%g) !>= 1.0\n", X[i * incX], scale);
      printf("%g * (%g) !>= 1.0\n", X[i * incX], 1.0/scale);
      printf("%g !>= 1.0\n", ratio);
      return 1;
    }
    if(ratio >= bound){
      printf("%g * (1.0/dscale(%g)) !< 2^diwidth()\n", X[i * incX], X[i * incX]);
      printf("%g * (1.0/%g) !< %g\n", X[i * incX], scale, bound);
      printf("%g * (%g) !< %g\n", X[i * incX], 1.0/scale, bound);
      printf("%g !< %g\n", ratio, bound);
      return 1;
    }
  }
  if(dindex(dscale(0.0)) != dindex(0.0)){
    printf("dindex(dscale(0.0)) != dindex(0.0)\n");
    return 1;
  }
  return 0;
}
