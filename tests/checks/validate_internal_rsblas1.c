#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>
#include "../common/test_vec.h"

static char* namebuf = "Validate rssum, rsasum, rsdot (1Big, 1BigPosNeg)";

const char* vecvec_name(int argc, char** argv) {
  return namebuf;
}

int vecvec_test(int argc, char** argv, int N, int incx, int incy) {
  int i;
	float small = 1.0 / 1024.0;			// 2^-10
	float big   = 1024.0 * 32;		// 2^15
  float ref;
  float refa;
  float res;

  vec_random_seed();

  //allocate vectors
  float *x    = svec_alloc(N, incx);
  float *y    = svec_alloc(N, incy);

  //fill empty space with random data to check increment
  svec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);
  svec_fill(N * incy, y, 1, vec_fill_RAND, 1.0, 1.0);

  //fill y with 1 where necessary
  svec_fill(N, y, incy, vec_fill_CONSTANT, 1.0, 1.0);

  //1 Big
  ref   = (N - 1) * small + big;

  //1 Big at beginning
  svec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[0]         = big;

  res = rssum(N, x, incx);
  if (res != ref) {
    printf("rssum(x) = %g != %g (1 Big at beginning)\n", res, ref);
    return 1;
  }

  res = rsasum(N, x, incx);
  if (res != ref) {
    printf("rsasum(x) = %g != %g (1 Big at beginning)\n", res, ref);
    return 1;
  }

  res = rsdot(N, x, incx, y, incy);
  if (res != ref) {
    printf("rsdot(x) = %g != %g (1 Big at beginning)\n", res, ref);
    return 1;
  }

  //1 Big in middle
  svec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[(N/2)*incx] = big;

  res = rssum(N, x, incx);
  if (res != ref) {
    printf("rssum(x) = %g != %g (1 Big in middle)\n", res, ref);
    return 1;
  }

  res = rsasum(N, x, incx);
  if (res != ref) {
    printf("rsasum(x) = %g != %g (1 Big in middle)\n", res, ref);
    return 1;
  }

  res = rsdot(N, x, incx, y, incy);
  if (res != ref) {
    printf("rsdot(x) = %g != %g (1 Big in middle)\n", res, ref);
    return 1;
  }

  //1 Big at end
  svec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[(N-1)*incx] = big;

  res = rssum(N, x, incx);
  if (res != ref) {
    printf("rssum(x) = %g != %g (1 Big at end)\n", res, ref);
    return 1;
  }

  res = rsasum(N, x, incx);
  if (res != ref) {
    printf("rsasum(x) = %g != %g (1 Big at end)\n", res, ref);
    return 1;
  }

  res = rsdot(N, x, incx, y, incy);
  if (res != ref) {
    printf("rsdot(x) = %g != %g (1 Big at end)\n", res, ref);
    return 1;
  }

  //1 Big pos neg
  ref   = (N - 2) * small;
  refa  = ((N - 2) * small) + 2 * big;

  //1 Big pos neg at beginning
  svec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[0]         = big;
  x[(N/2)*incx] = -big;

  res = rssum(N, x, incx);
  if (res != ref) {
    printf("rssum(x) = %g != %g (1 Big pos neg at beginning)\n", res, ref);
    return 1;
  }

  res = rsasum(N, x, incx);
  if (res != refa) {
    printf("rsasum(x) = %g != %g (1 Big pos neg at beginning)\n", res, refa);
    return 1;
  }

  res = rsdot(N, x, incx, y, incy);
  if (res != ref) {
    printf("rsdot(x) = %g != %g (1 Big pos neg at beginning)\n", res, ref);
    return 1;
  }

  //1 Big pos neg at ends
  svec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[0]         = big;
  x[(N-1)*incx] = -big;

  res = rssum(N, x, incx);
  if (res != ref) {
    printf("rssum(x) = %g != %g (1 Big pos neg at ends)\n", res, ref);
    return 1;
  }

  res = rsasum(N, x, incx);
  if (res != refa) {
    printf("rsasum(x) = %g != %g (1 Big pos neg at ends)\n", res, refa);
    return 1;
  }

  res = rsdot(N, x, incx, y, incy);
  if (res != ref) {
    printf("rsdot(x) = %g != %g (1 Big pos neg at ends)\n", res, ref);
    return 1;
  }

  //1 Big pos neg at end
  svec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[(N/2)*incx] = big;
  x[(N-1)*incx] = -big;

  res = rssum(N, x, incx);
  if (res != ref) {
    printf("rssum(x) = %g != %g (1 Big pos neg at end)\n", res, ref);
    return 1;
  }

  res = rsasum(N, x, incx);
  if (res != refa) {
    printf("rsasum(x) = %g != %g (1 Big pos neg at end)\n", res, refa);
    return 1;
  }

  res = rsdot(N, x, incx, y, incy);
  if (res != ref) {
    printf("rsdot(x) = %g != %g (1 Big pos neg at end)\n", res, ref);
    return 1;
  }

  free(x);
  free(y);

  return 0;
}
