#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>
#include "../common/test_vec.h"

static char* namebuf = "Validate rcsum, rscasum, rcdotu/c (1Big, 1BigPosNeg)";

const char* vecvec_name(int argc, char** argv) {
  return namebuf;
}

int vecvec_check(int argc, char** argv, int N, int incx, int incy) {
  int i;
  float small = 1.0 / 1024.0;			// 2^-10
  float big   = 1024.0 * 32;		// 2^15
  float complex ref;
  float complex refa;
  float complex res;

  vec_random_seed();

  //allocate vectors
  float complex *x    = cvec_alloc(N, incx);
  float complex *y    = cvec_alloc(N, incx);

  //fill empty space with random data to check increment
  cvec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);
  cvec_fill(N * incy, y, 1, vec_fill_RAND, 1.0, 1.0);

  //fill y with -i where necessary
  cvec_fill(N, y, incy, vec_fill_CONSTANT, -_Complex_I, 1.0);

  //1 Big
  ref   = (N - 1) * small + big;

  //1 Big at beginning
  cvec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[0]         = big;

  res = rcsum(N, x, incx);
  if (res != ref) {
    printf("rcsum(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rscasum(N, x, incx);
  if (res != ref) {
    printf("rscasum(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rcdotu(N, x, incx, y, incy);
  if (res != ref * -_Complex_I) {
    printf("rcdotu(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * -_Complex_I), CIMAG_(ref * -_Complex_I));
    return 1;
  }

  res = rcdotc(N, x, incx, y, incy);
  if (res != ref * _Complex_I) {
    printf("rcdotc(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * _Complex_I), CIMAG_(ref * _Complex_I));
    return 1;
  }

  //1 Big in middle
  cvec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[(N/2)*incx] = big;

  res = rcsum(N, x, incx);
  if (res != ref) {
    printf("rcsum(x) = %g + %gi != %g + %gi (1 Big in middle)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rscasum(N, x, incx);
  if (res != ref) {
    printf("rscasum(x) = %g + %gi != %g + %gi (1 Big in middle)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rcdotu(N, x, incx, y, incy);
  if (res != ref * -_Complex_I) {
    printf("rcdotu(x) = %g + %gi != %g + %gi (1 Big in middle)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rcdotc(N, x, incx, y, incy);
  if (res != ref * _Complex_I) {
    printf("rcdotc(x) = %g + %gi != %g + %gi (1 Big in middle)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  //1 Big at end
  cvec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[(N-1)*incx] = big;

  res = rcsum(N, x, incx);
  if (res != ref) {
    printf("rcsum(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rscasum(N, x, incx);
  if (res != ref) {
    printf("rscasum(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rcdotu(N, x, incx, y, incy);
  if (res != ref * -_Complex_I) {
    printf("rcdotu(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * -_Complex_I), CIMAG_(ref * -_Complex_I));
    return 1;
  }

  res = rcdotc(N, x, incx, y, incy);
  if (res != ref * _Complex_I) {
    printf("rcdotc(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * _Complex_I), CIMAG_(ref * _Complex_I));
    return 1;
  }

  //1 Big pos neg
  ref   = (N - 2) * small;
  refa  = ((N - 2) * small) + 2 * big;

  //1 Big pos neg at beginning
  cvec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[0]         = big;
  x[(N/2)*incx] = -big;

  res = rcsum(N, x, incx);
  if (res != ref) {
    printf("rcsum(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rscasum(N, x, incx);
  if (res != refa) {
    printf("rscasum(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(refa));
    return 1;
  }

  res = rcdotu(N, x, incx, y, incy);
  if (res != ref * -_Complex_I) {
    printf("rcdotu(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * -_Complex_I), CIMAG_(ref * -_Complex_I));
    return 1;
  }

  res = rcdotc(N, x, incx, y, incy);
  if (res != ref * _Complex_I) {
    printf("rcdotc(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * _Complex_I), CIMAG_(ref * _Complex_I));
    return 1;
  }

  //1 Big pos neg at ends
  cvec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[0]         = big;
  x[(N-1)*incx] = -big;

  res = rcsum(N, x, incx);
  if (res != ref) {
    printf("rcsum(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rscasum(N, x, incx);
  if (res != refa) {
    printf("rscasum(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(refa));
    return 1;
  }

  res = rcdotu(N, x, incx, y, incy);
  if (res != ref * -_Complex_I) {
    printf("rcdotu(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * -_Complex_I), CIMAG_(ref * -_Complex_I));
    return 1;
  }

  res = rcdotc(N, x, incx, y, incy);
  if (res != ref * _Complex_I) {
    printf("rcdotc(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * _Complex_I), CIMAG_(ref * _Complex_I));
    return 1;
  }

  //1 Big pos neg at end
  cvec_fill(N, x, incx, vec_fill_CONSTANT, small, 1.0);
  x[(N/2)*incx] = big;
  x[(N-1)*incx] = -big;

  res = rcsum(N, x, incx);
  if (res != ref) {
    printf("rcsum(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rscasum(N, x, incx);
  if (res != refa) {
    printf("rscasum(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(refa));
    return 1;
  }

  res = rcdotu(N, x, incx, y, incy);
  if (res != ref * -_Complex_I) {
    printf("rcdotu(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * -_Complex_I), CIMAG_(ref * -_Complex_I));
    return 1;
  }

  res = rcdotc(N, x, incx, y, incy);
  if (res != ref * _Complex_I) {
    printf("rcdotc(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref * _Complex_I), CIMAG_(ref * _Complex_I));
    return 1;
  }

  free(x);
  free(y);

  return 0;
}
