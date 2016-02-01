#include <idxdBLAS.h>
#include <idxd.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/test_vecvec_header.h"

int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  (void)argc;
  (void)argv;
  return "Validate idxdBLAS_damax(m) X=+big";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  double small = 1.0 / (1024.0 * 1024.0 * 128.0); // 2^-27
  double big   = 1024.0 * 1024.0 * 128.0;         // 2^27
  double ref;
  double res;

  util_random_seed();

  //allocate vectors
  double *X    = util_dvec_alloc(N, incX);
  double *Y    = util_dvec_alloc(N, incY);

  //fill Y with 1 where necessary
  util_dvec_fill(N, Y, incY, util_Vec_Constant, 1, 1.0);

  //1 Big
  ref   = big;

  //1 Big at beginning
  util_dvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = -big;

  res = idxdBLAS_damax(N, X, incX);
  if (res != ref) {
    printf("idxdBLAS_damax(X) = %g != %g (1 Big at beginning)\n", res, ref);
    return 1;
  }

  res = idxdBLAS_damaxm(N, X, incX, Y, incY);
  if (res != ref) {
    printf("idxdBLAS_damaxm(X) = %g != %g (1 Big at beginning)\n", res, ref);
    return 1;
  }

  //1 Big at end
  util_dvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N-1)*incX]         = -big;

  res = idxdBLAS_damax(N, X, incX);
  if (res != ref) {
    printf("idxdBLAS_damax(X) = %g != %g (1 Big at end)\n", res, ref);
    return 1;
  }

  res = idxdBLAS_damaxm(N, X, incX, Y, incY);
  if (res != ref) {
    printf("idxdBLAS_damaxm(X) = %g != %g (1 Big at end)\n", res, ref);
    return 1;
  }

  free(X);
  free(Y);

  return 0;
}
