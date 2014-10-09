#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>
#include "../common/test_vec.h"

static char* namebuf = "Validate rzsum, rdzasum, rzdotu/c (1Big, 1BigPosNeg)";

const char* vecvec_name(int argc, char** argv) {
  return namebuf;
}

int vecvec_check(int argc, char** argv, int N, int incx, int incy) {
  int i;
  double small = 1.0 / (1024.0 * 1024.0);       // 2^-20
  double big   = 1024.0 * 1024.0 * 1024.0 * 32; // 2^35
  double complex ref;
  double complex refa;
  double complex res;

  vec_random_seed();

  //allocate vectors
  double complex *x    = zvec_alloc(N, incx);
  double complex *y    = zvec_alloc(N, incy);

  //fill empty space with random data to check increment
  zvec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);
  zvec_fill(N * incy, y, 1, vec_fill_RAND, 1.0, 1.0);

  //fill y with -i where necessary
  zvec_fill(N, y, incy, vec_fill_CONSTANT, -_Complex_I, 1.0);

  //1 Big
  ref   = (N - 1) * small + big;

  //1 Big at beginning
  zvec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[0]         = big;

  res = rzsum(N, x, incx);
  if (res != ref) {
    printf("rzsum(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rdzasum(N, x, incx);
  if (res != ref) {
    printf("rdzasum(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rzdotu(N, x, incx, y, incy);
  if (res != ref * -_Complex_I) {
    printf("rzdotu(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * -_Complex_I), CIMAG_(ref * -_Complex_I));
    return 1;
  }

  res = rzdotc(N, x, incx, y, incy);
  if (res != ref * _Complex_I) {
    printf("rzdotc(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * _Complex_I), CIMAG_(ref * _Complex_I));
    return 1;
  }

  //1 Big in middle
  zvec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[(N/2)*incx] = big;

  res = rzsum(N, x, incx);
  if (res != ref) {
    printf("rzsum(x) = %g + %gi != %g + %gi (1 Big in middle)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rdzasum(N, x, incx);
  if (res != ref) {
    printf("rdzasum(x) = %g + %gi != %g + %gi (1 Big in middle)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rzdotu(N, x, incx, y, incy);
  if (res != ref * -_Complex_I) {
    printf("rzdotu(x) = %g + %gi != %g + %gi (1 Big in middle)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * -_Complex_I), CIMAG_(ref * -_Complex_I));
    return 1;
  }

  res = rzdotc(N, x, incx, y, incy);
  if (res != ref * _Complex_I) {
    printf("rzdotc(x) = %g + %gi != %g + %gi (1 Big in middle)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * _Complex_I), CIMAG_(ref * _Complex_I));
    return 1;
  }

  //1 Big at end
  zvec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[(N-1)*incx] = big;

  res = rzsum(N, x, incx);
  if (res != ref) {
    printf("rzsum(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rdzasum(N, x, incx);
  if (res != ref) {
    printf("rdzasum(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rzdotu(N, x, incx, y, incy);
  if (res != ref * -_Complex_I) {
    printf("rzdotu(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * -_Complex_I), CIMAG_(ref * -_Complex_I));
    return 1;
  }

  res = rzdotc(N, x, incx, y, incy);
  if (res != ref * _Complex_I) {
    printf("rzdotc(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * _Complex_I), CIMAG_(ref * _Complex_I));
    return 1;
  }

  //1 Big pos neg
  ref   = (N - 2) * small;
  refa  = ((N - 2) * small) + 2 * big;

  //1 Big pos neg at beginning
  zvec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[0]         = big;
  x[(N/2)*incx] = -big;

  res = rzsum(N, x, incx);
  if (res != ref) {
    printf("rzsum(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rdzasum(N, x, incx);
  if (res != refa) {
    printf("rdzasum(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(refa));
    return 1;
  }

  res = rzdotu(N, x, incx, y, incy);
  if (res != ref * -_Complex_I) {
    printf("rzdotu(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * -_Complex_I), CIMAG_(ref * -_Complex_I));
    return 1;
  }

  res = rzdotc(N, x, incx, y, incy);
  if (res != ref * _Complex_I) {
    printf("rzdotc(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * _Complex_I), CIMAG_(ref * _Complex_I));
    return 1;
  }

  //1 Big pos neg at ends
  zvec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[0]         = big;
  x[(N-1)*incx] = -big;

  res = rzsum(N, x, incx);
  if (res != ref) {
    printf("rzsum(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rdzasum(N, x, incx);
  if (res != refa) {
    printf("rdzasum(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(refa));
    return 1;
  }

  res = rzdotu(N, x, incx, y, incy);
  if (res != ref * -_Complex_I) {
    printf("rzdotu(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * -_Complex_I), CIMAG_(ref * -_Complex_I));
    return 1;
  }

  res = rzdotc(N, x, incx, y, incy);
  if (res != ref * _Complex_I) {
    printf("rzdotc(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * _Complex_I), CIMAG_(ref * _Complex_I));
    return 1;
  }

  //1 Big pos neg at end
  zvec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[(N/2)*incx] = big;
  x[(N-1)*incx] = -big;

  res = rzsum(N, x, incx);
  if (res != ref) {
    printf("rzsum(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rdzasum(N, x, incx);
  if (res != refa) {
    printf("rdzasum(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(refa));
    return 1;
  }

  res = rzdotu(N, x, incx, y, incy);
  if (res != ref * -_Complex_I) {
    printf("rzdotu(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * -_Complex_I), CIMAG_(ref * -_Complex_I));
    return 1;
  }

  res = rzdotc(N, x, incx, y, incy);
  if (res != ref * _Complex_I) {
    printf("rzdotc(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * _Complex_I), CIMAG_(ref * _Complex_I));
    return 1;
  }

  free(x);
  free(y);

  return 0;
}
