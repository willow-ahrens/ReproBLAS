#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/test_vecvec_header.h"

int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  return "Validate rcsum, rscasum, rcdotu/c (1Big, 1BigPosNeg)";
}

int vecvec_test(int argc, char** argv, int N, int incx, int incy) {
  int i;
  float small = 1.0 / 1024.0; // 2^-10
  float big   = 1024.0 * 32;  // 2^15
  float complex ref;
  float complex refa;
  float complex refdu;
  float complex refdc;
  float complex res;

  util_random_seed();

  //allocate vectors
  float complex *x    = util_cvec_alloc(N, incx);
  float complex *y    = util_cvec_alloc(N, incx);

  //fill y with -i where necessary
  util_cvec_fill(N, y, incy, util_Vec_Constant, -_Complex_I, 1.0);

  //1 Big
  ref   = (N - 1) * small + big;
  refa  = ref;
  refdu  = ref * -_Complex_I;
  refdc  = ref * _Complex_I;

  //1 Big at beginning
  util_cvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[0]         = big;

  res = rcsum(N, x, incx);
  if (res != ref) {
    printf("rcsum(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rscasum(N, x, incx);
  if (res != refa) {
    printf("rscasum(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(refa), CIMAG_(refa));
    return 1;
  }

  res = rcdotu(N, x, incx, y, incy);
  if (res != refdu) {
    printf("rcdotu(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(refdu), CIMAG_(refdu));
    return 1;
  }

  res = rcdotc(N, x, incx, y, incy);
  if (res != refdc) {
    printf("rcdotc(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(refdc), CIMAG_(refdc));
    return 1;
  }

  //1 Big in middle
  util_cvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[(N/2)*incx] = big;

  res = rcsum(N, x, incx);
  if (res != ref) {
    printf("rcsum(x) = %g + %gi != %g + %gi (1 Big in middle)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rscasum(N, x, incx);
  if (res != refa) {
    printf("rscasum(x) = %g + %gi != %g + %gi (1 Big in middle)\n", CREAL_(res), CIMAG_(res), CREAL_(refa), CIMAG_(refa));
    return 1;
  }

  res = rcdotu(N, x, incx, y, incy);
  if (res != refdu) {
    printf("rcdotu(x) = %g + %gi != %g + %gi (1 Big in middle)\n", CREAL_(res), CIMAG_(res), CREAL_(refdu), CIMAG_(refdu));
    return 1;
  }

  res = rcdotc(N, x, incx, y, incy);
  if (res != refdc) {
    printf("rcdotc(x) = %g + %gi != %g + %gi (1 Big in middle)\n", CREAL_(res), CIMAG_(res), CREAL_(refdc), CIMAG_(refdc));
    return 1;
  }

  //1 Big at end
  util_cvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[(N-1)*incx] = big;

  res = rcsum(N, x, incx);
  if (res != ref) {
    printf("rcsum(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rscasum(N, x, incx);
  if (res != refa) {
    printf("rscasum(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(refa), CIMAG_(refa));
    return 1;
  }

  res = rcdotu(N, x, incx, y, incy);
  if (res != refdu) {
    printf("rcdotu(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(refdu), CIMAG_(refdu));
    return 1;
  }

  res = rcdotc(N, x, incx, y, incy);
  if (res != refdc) {
    printf("rcdotc(x) = %g + %gi != %g + %gi (1 Big at end)\n", CREAL_(res), CIMAG_(res), CREAL_(refdc), CIMAG_(refdc));
    return 1;
  }

  //1 Big pos neg
  ref   = (N - 2) * small;
  refa  = ((N - 2) * small) + 2 * big;
  refdu  = ref * -_Complex_I;
  refdc  = ref * _Complex_I;

  //1 Big pos neg at beginning
  util_cvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[0]         = big;
  x[(N/2)*incx] = -big;

  res = rcsum(N, x, incx);
  if (res != ref) {
    printf("rcsum(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rscasum(N, x, incx);
  if (res != refa) {
    printf("rscasum(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(refa), CIMAG_(refa));
    return 1;
  }

  res = rcdotu(N, x, incx, y, incy);
  if (res != refdu) {
    printf("rcdotu(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(refdu), CIMAG_(refdu));
    return 1;
  }

  res = rcdotc(N, x, incx, y, incy);
  if (res != refdc) {
    printf("rcdotc(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", CREAL_(res), CIMAG_(res), CREAL_(refdc), CIMAG_(refdc));
    return 1;
  }

  //1 Big pos neg at ends
  util_cvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[0]         = big;
  x[(N-1)*incx] = -big;

  res = rcsum(N, x, incx);
  if (res != ref) {
    printf("rcsum(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rscasum(N, x, incx);
  if (res != refa) {
    printf("rscasum(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", CREAL_(res), CIMAG_(res), CREAL_(refa), CIMAG_(refa));
    return 1;
  }

  res = rcdotu(N, x, incx, y, incy);
  if (res != refdu) {
    printf("rcdotu(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", CREAL_(res), CIMAG_(res), CREAL_(refdu), CIMAG_(refdu));
    return 1;
  }

  res = rcdotc(N, x, incx, y, incy);
  if (res != refdc) {
    printf("rcdotc(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", CREAL_(res), CIMAG_(res), CREAL_(refdc), CIMAG_(refdc));
    return 1;
  }

  //1 Big pos neg at end
  util_cvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[(N/2)*incx] = big;
  x[(N-1)*incx] = -big;

  res = rcsum(N, x, incx);
  if (res != ref) {
    printf("rcsum(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
    return 1;
  }

  res = rscasum(N, x, incx);
  if (res != refa) {
    printf("rscasum(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", CREAL_(res), CIMAG_(res), CREAL_(refa), CIMAG_(refa));
    return 1;
  }

  res = rcdotu(N, x, incx, y, incy);
  if (res != refdu) {
    printf("rcdotu(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", CREAL_(res), CIMAG_(res), CREAL_(refdu), CIMAG_(refdu));
    return 1;
  }

  res = rcdotc(N, x, incx, y, incy);
  if (res != refdc) {
    printf("rcdotc(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", CREAL_(res), CIMAG_(res), CREAL_(refdc), CIMAG_(refdc));
    return 1;
  }

  free(x);
  free(y);

  return 0;
}
