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
  return "Validate samax(m) (1Big)";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  float small = 1.0 / 1024.0;			// 2^-10
  float big   = 1024.0 * 32;		// 2^15
  float ref;
  float res;

  util_random_seed();

  //allocate vectors
  float *X    = util_svec_alloc(N, incX);
  float *Y    = util_svec_alloc(N, incY);

  //fill Y with 1 where necessary
  util_svec_fill(N, Y, incY, util_Vec_Constant, 1, 1.0);

  //1 Big
  ref   = big;

  //1 Big at beginning
  util_svec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = -big;

  res = samax(N, X, incX);
  if (res != ref) {
    printf("samax(X) = %g != %g (1 Big at beginning)\n", res, ref);
    return 1;
  }

  res = samaxm(N, X, incX, Y, incY);
  if (res != ref) {
    printf("samaxm(X) = %g != %g (1 Big at beginning)\n", res, ref);
    return 1;
  }

  //1 Big at end
  util_svec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N-1)*incX]         = -big;

  res = samax(N, X, incX);
  if (res != ref) {
    printf("samax(X) = %g != %g (1 Big at end)\n", res, ref);
    return 1;
  }

  res = samaxm(N, X, incX, Y, incY);
  if (res != ref) {
    printf("samaxm(X) = %g != %g (1 Big at end)\n", res, ref);
    return 1;
  }

  free(X);
  free(Y);

  return 0;
}
