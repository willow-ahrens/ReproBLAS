#include <idxdBLAS.h>
#include <idxd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../common/test_vecvec_header.h"
#include "../common/test_util.h"


int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  (void)argc;
  (void)argv;
  return "Verify dmindex";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  (void)incY;
  int i;
  int index;

  util_random_seed();

  //allocate vector
  double *X = util_dvec_alloc(N * ((DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH) + 1, incX);

  //check
  for (i = 0; i < N * ((DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH) + 1; i++) {
    X[i * incX] = *idxd_dmbins(i/N) * (1.25/1.5 + 0.5/1.5 * util_drand48());
  }
  for (i = 0; i < N * ((DBL_MAX_EXP - DBL_MIN_EXP)/DIWIDTH) + 1; i++) {
    index = idxd_dmindex(X + i * incX);
    if (*idxd_dmbins(index)*(1.25/1.5) > X[i * incX]){
      printf("X < 1.25 * 2^(i * W - 1)\n");
      printf("%g < 1.25 * %g\n", X[i * incX], *idxd_dmbins(index));
      return 1;
    }
    if (*idxd_dmbins(index)*(1.75/1.5) < X[i * incX]){
      printf("X > 1.75 * 2^(i * W - 1)\n");
      printf("%g > 1.75 * %g\n", X[i * incX], *idxd_dmbins(index));
      return 1;
    }
  }
  free(X);
  return 0;
}
