#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/test_vecvec_header.h"

int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  return "Validate zamax(m) (1Big)";
}

int vecvec_test(int argc, char** argv, int N, int incx, int incy) {
  int i;
  double small = 1.0 / (1024.0 * 1024.0);       // 2^-20
  double big   = 1024.0 * 1024.0 * 1024.0 * 32; // 2^35
  double complex ref;
  double complex res;

  util_random_seed();

  //allocate vectors
  double complex *x    = util_zvec_alloc(N, incx);
  double complex *y    = util_zvec_alloc(N, incy);

  //fill y with 1 where necessary
  util_zvec_fill(N, y, incy, util_Vec_Constant, 1.0, 1.0);

  //1 Big
  ref   = big + _Complex_I * big;

  //1 Big at beginning
  util_zvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[0]         = -big + -_Complex_I * big;

  res = zamax(N, x, incx);
  if (res != ref) {
    printf("zamax(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  res = zamaxm(N, x, incx, y, incy);
  if (res != ref) {
    printf("zamaxm(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  //1 Big at end
  util_zvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[(N-1)*incx]         = -big + -_Complex_I * big;

  res = zamax(N, x, incx);
  if (res != ref) {
    printf("zamax(x) = %g + %gi != %g + %gi (1 Big at end)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  res = zamaxm(N, x, incx, y, incy);
  if (res != ref) {
    printf("zamaxm(x) = %g + %gi != %g + %gi (1 Big at end)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  free(x);
  free(y);

  return 0;
}
