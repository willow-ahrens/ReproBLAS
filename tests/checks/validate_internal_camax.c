#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/test_vecvec_header.h"

int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  return "Validate camax(m) (1Big)";
}

int vecvec_test(int argc, char** argv, int N, int incx, int incy) {
  int i;
  float small = 1.0 / 1024.0;			// 2^-10
  float big   = 1024.0 * 32;		// 2^15
  float complex ref;
  float complex res;

  util_random_seed();

  //allocate vectors
  float complex *x    = util_cvec_alloc(N, incx);
  float complex *y    = util_cvec_alloc(N, incy);

  //fill y with 1 where necessary
  util_cvec_fill(N, y, incy, util_Vec_Constant, 1.0, 1.0);

  //1 Big
  ref   = big + _Complex_I * big;

  //1 Big at beginning
  util_cvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[0]         = -big + -_Complex_I * big;

  res = camax(N, x, incx);
  if (res != ref) {
    printf("camax(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = camaxm(N, x, incx, y, incy);
  if (res != ref) {
    printf("camaxm(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  //1 Big at end
  util_cvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[(N-1)*incx]         = -big + -_Complex_I * big;

  res = camax(N, x, incx);
  if (res != ref) {
    printf("camax(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = camaxm(N, x, incx, y, incy);
  if (res != ref) {
    printf("camaxm(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  free(x);
  free(y);

  return 0;
}
