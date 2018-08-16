#include <binnedBLAS.h>
#include <binned.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/test_vecvec_header.h"

int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  (void)argc;
  (void)argv;
  return "Validate camax(m) X=+big";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  float small = 1.0 / (1024.0 * 4.0); // 2^-12
  float big   = 1024.0 * 8.0;         // 2^13
  float complex ref;
  float complex res;

  util_random_seed();

  //allocate vectors
  float complex *X    = util_cvec_alloc(N, incX);
  float complex *Y    = util_cvec_alloc(N, incY);

  //fill Y with 1 where necessary
  util_cvec_fill(N, Y, incY, util_Vec_Constant, 1.0, 1.0);

  //1 Big
  ref   = big + _Complex_I * big;

  //1 Big at beginning
  util_cvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = -big + -_Complex_I * big;

  binnedBLAS_camax_sub(N, X, incX, &res);
  if (res != ref) {
    printf("camax(X) = %g + %gi != %g + %gi (1 Big at beginning)\n", crealf(res), cimagf(res), crealf(ref), cimagf(ref));
    return 1;
  }

  binnedBLAS_camaxm_sub(N, X, incX, Y, incY, &res);
  if (res != ref) {
    printf("camaxm(X) = %g + %gi != %g + %gi (1 Big at beginning)\n", crealf(res), cimagf(res), crealf(ref), cimagf(ref));
    return 1;
  }

  //1 Big at end
  util_cvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N-1)*incX]         = -big + -_Complex_I * big;

  binnedBLAS_camax_sub(N, X, incX, &res);
  if (res != ref) {
    printf("camax(X) = %g + %gi != %g + %gi (1 Big at end)\n", crealf(res), cimagf(res), crealf(ref), cimagf(ref));
    return 1;
  }

  binnedBLAS_camaxm_sub(N, X, incX, Y, incY, &res);
  if (res != ref) {
    printf("camaxm(X) = %g + %gi != %g + %gi (1 Big at end)\n", crealf(res), cimagf(res), crealf(ref), cimagf(ref));
    return 1;
  }

  free(X);
  free(Y);

  return 0;
}
