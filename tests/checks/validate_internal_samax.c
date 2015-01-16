#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/test_vecvec_header.h"

int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  return "Validate samax(m) (1Big)";
}

int vecvec_test(int argc, char** argv, int N, int incx, int incy) {
  int i;
  float small = 1.0 / 1024.0;			// 2^-10
  float big   = 1024.0 * 32;		// 2^15
  float ref;
  float res;

  vec_random_seed();

  //allocate vectors
  float *x    = svec_alloc(N, incx);
  float *y    = svec_alloc(N, incy);

  //fill empty space with random data to check increment
  svec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);
  svec_fill(N * incy, y, 1, vec_fill_RAND, 1.0, 1.0);

  //fill y with 1 where necessary
  svec_fill(N, y, incy, vec_fill_CONSTANT, 1, 1.0);

  //1 Big
  ref   = big;

  //1 Big at beginning
  svec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[0]         = -big;

  res = samax(N, x, incx);
  if (res != ref) {
    printf("samax(x) = %g != %g (1 Big at beginning)\n", res, ref);
    return 1;
  }

  res = samaxm(N, x, incx, y, incy);
  if (res != ref) {
    printf("samaxm(x) = %g != %g (1 Big at beginning)\n", res, ref);
    return 1;
  }

  //1 Big at end
  svec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[(N-1)*incx]         = -big;

  res = samax(N, x, incx);
  if (res != ref) {
    printf("samax(x) = %g != %g (1 Big at end)\n", res, ref);
    return 1;
  }

  res = samaxm(N, x, incx, y, incy);
  if (res != ref) {
    printf("samaxm(x) = %g != %g (1 Big at end)\n", res, ref);
    return 1;
  }

  free(x);
  free(y);

  return 0;
}
