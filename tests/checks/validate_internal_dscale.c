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

  util_random_seed();

  //allocate vector
  double *X = util_dvec_alloc(N*(DBL_MAX_EXP - DBL_MIN_EXP), incX);

  //check
  for (i = 0; i < (DBL_MAX_EXP - DBL_MIN_EXP) * N; i++) {
    X[i * incX] = ldexp(0.5 + 0.5 * util_drand48(), (i/N) + DBL_MIN_EXP);
  }
  for (i = 0; i < N * (DBL_MAX_EXP - DBL_MIN_EXP); i++) {
    if (fabs(dindex(X[i * incX]/dscale(X[i * incX])) - dindex(1.0)) > 1){
      printf("|dindex(%g/dscale(%g)) - dindex(1.0)| > 1\n", X[i * incX], X[i * incX]);
      printf("|dindex(%g/%g) - dindex(1.0)| > 1\n", X[i * incX], dscale(X[i * incX]));
      printf("|dindex(%g) - dindex(1.0)| > 1\n", X[i * incX]/dscale(X[i * incX]));
      printf("|%d - %d| > 1\n", dindex(X[i * incX]/dscale(X[i * incX])), dindex(1.0));
      return 1;
    }
  }
  if(dindex(dscale(0.0)) != dindex(0.0)){
    printf("dindex(dscale(0.0)) != dindex(0.0)\n");
    return 1;
  }
  return 0;
}
