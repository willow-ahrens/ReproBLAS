#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/test_vecvec_header.h"

int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  (void)argc;
  (void)argv;
  return "Validate ufp";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  (void)incY;
  int i;
  double bin = 1.0;
  double ufpX;

  util_random_seed();

  //allocate vector
  double *X = util_dvec_alloc(N, incX);

  //check
  for (i = 0; i < N; i++) {
    X[i * incX] = 3 * bin;
    bin *= 2;

    ufpX = ufp(X[i * incX]);
    if (ufpX != bin) {
      printf("ufp(%g) = %g != %g\n", X[i * incX], ufpX, bin);
      return 1;
    }
  }
  return 0;
}
