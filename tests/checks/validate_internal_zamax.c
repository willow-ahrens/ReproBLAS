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
  return "Validate zamax(m) (1Big)";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  double small = 1.0 / (1024.0 * 1024.0);       // 2^-20
  double big   = 1024.0 * 1024.0 * 1024.0 * 32; // 2^35
  double complex ref;
  double complex res;

  util_random_seed();

  //allocate vectors
  double complex *X    = util_zvec_alloc(N, incX);
  double complex *Y    = util_zvec_alloc(N, incY);

  //fill Y with 1 where necessary
  util_zvec_fill(N, Y, incY, util_Vec_Constant, 1.0, 1.0);

  //1 Big
  ref   = big + _Complex_I * big;

  //1 Big at beginning
  util_zvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = -big + -_Complex_I * big;

  zamax_sub(N, X, incX, &res);
  if (res != ref) {
    printf("zamax(X) = %g + %gi != %g + %gi (1 Big at beginning)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  zamaxm_sub(N, X, incX, Y, incY, &res);
  if (res != ref) {
    printf("zamaxm(X) = %g + %gi != %g + %gi (1 Big at beginning)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  //1 Big at end
  util_zvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N-1)*incX]         = -big + -_Complex_I * big;

  zamax_sub(N, X, incX, &res);
  if (res != ref) {
    printf("zamax(X) = %g + %gi != %g + %gi (1 Big at end)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  zamaxm_sub(N, X, incX, Y, incY, &res);
  if (res != ref) {
    printf("zamaxm(X) = %g + %gi != %g + %gi (1 Big at end)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  free(X);
  free(Y);

  return 0;
}
