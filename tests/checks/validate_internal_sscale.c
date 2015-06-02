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
  return "Validate sscale internally";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  (void)incY;
  int i;

  util_random_seed();

  //allocate vector
  float *X = util_svec_alloc(N*(FLT_MAX_EXP - FLT_MIN_EXP), incX);

  //check
  for (i = 0; i < (FLT_MAX_EXP - FLT_MIN_EXP) * N; i++) {
    X[i * incX] = ldexpf(0.5 + 0.5 * util_drand48(), (i/N) + FLT_MIN_EXP);
  }
  for (i = 0; i < N * (FLT_MAX_EXP - FLT_MIN_EXP); i++) {
    if (abs(sindex(X[i * incX]/sscale(X[i * incX])) - sindex(1.0)) > 1){
      printf("|sindex(%g/sscale(%g)) - sindex(1.0)| > 1\n", X[i * incX], X[i * incX]);
      printf("|sindex(%g/%g) - sindex(1.0)| > 1\n", X[i * incX], sscale(X[i * incX]));
      printf("|sindex(%g) - sindex(1.0)| > 1\n", X[i * incX]/sscale(X[i * incX]));
      printf("|%d - %d| > 1\n", sindex(X[i * incX]/sscale(X[i * incX])), sindex(1.0));
      return 1;
    }
  }
  if(sindex(sscale(0.0)) != sindex(0.0)){
    printf("sindex(sscale(0.0)) != sindex(0.0)\n");
    return 1;
  }
  return 0;
}
