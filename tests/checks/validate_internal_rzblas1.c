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
  return "Validate rzsum, rdzasum, rzdotu/c (1Big, 1BigPosNeg)";
}

int vecvec_test(int argc, char** argv, int N, int incx, int incy) {
  (void)argc;
  (void)argv;
  double small = 1.0 / (1024.0 * 1024.0);       // 2^-20
  double big   = 1024.0 * 1024.0 * 1024.0 * 32; // 2^35
  double complex ref;
  double complex refa;
  double complex refdu;
  double complex refdc;
  double complex res;

  util_random_seed();

  //allocate vectors
  double complex *x    = util_zvec_alloc(N, incx);
  double complex *y    = util_zvec_alloc(N, incy);

  //fill y with -i where necessary
  util_zvec_fill(N, y, incy, util_Vec_Constant, -_Complex_I, 1.0);

  //1 Big
  ref   = (N - 1) * small + big;
  refa  = ref;
  refdu  = ref * -_Complex_I;
  refdc  = ref * _Complex_I;

  //1 Big at beginning
  util_zvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[0]         = big;

  res = rzsum(N, x, incx);
  if (res != ref) {
    printf("rzsum(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  res = rdzasum(N, x, incx);
  if (res != refa) {
    printf("rdzasum(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refa), ZIMAG_(refa));
    return 1;
  }

  res = rzdotu(N, x, incx, y, incy);
  if (res != refdu) {
    printf("rzdotu(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refdu), ZIMAG_(refdu));
    return 1;
  }

  res = rzdotc(N, x, incx, y, incy);
  if (res != refdc) {
    printf("rzdotc(x) = %g + %gi != %g + %gi (1 Big at beginning)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refdc), ZIMAG_(refdc));
    return 1;
  }

  //1 Big in middle
  util_zvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[(N/2)*incx] = big;

  res = rzsum(N, x, incx);
  if (res != ref) {
    printf("rzsum(x) = %g + %gi != %g + %gi (1 Big in middle)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  res = rdzasum(N, x, incx);
  if (res != refa) {
    printf("rdzasum(x) = %g + %gi != %g + %gi (1 Big in middle)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refa), ZIMAG_(refa));
    return 1;
  }

  res = rzdotu(N, x, incx, y, incy);
  if (res != refdu) {
    printf("rzdotu(x) = %g + %gi != %g + %gi (1 Big in middle)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refdu), ZIMAG_(refdu));
    return 1;
  }

  res = rzdotc(N, x, incx, y, incy);
  if (res != refdc) {
    printf("rzdotc(x) = %g + %gi != %g + %gi (1 Big in middle)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refdc), ZIMAG_(refdc));
    return 1;
  }

  //1 Big at end
  util_zvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[(N-1)*incx] = big;

  res = rzsum(N, x, incx);
  if (res != ref) {
    printf("rzsum(x) = %g + %gi != %g + %gi (1 Big at end)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  res = rdzasum(N, x, incx);
  if (res != refa) {
    printf("rdzasum(x) = %g + %gi != %g + %gi (1 Big at end)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refa), ZIMAG_(refa));
    return 1;
  }

  res = rzdotu(N, x, incx, y, incy);
  if (res != refdu) {
    printf("rzdotu(x) = %g + %gi != %g + %gi (1 Big at end)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refdu), ZIMAG_(refdu));
    return 1;
  }

  res = rzdotc(N, x, incx, y, incy);
  if (res != refdc) {
    printf("rzdotc(x) = %g + %gi != %g + %gi (1 Big at end)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refdc), ZIMAG_(refdc));
    return 1;
  }

  //1 Big pos neg
  ref   = (N - 2) * small;
  refa  = ((N - 2) * small) + 2 * big;
  refdu  = ref * -_Complex_I;
  refdc  = ref * _Complex_I;

  //1 Big pos neg at beginning
  util_zvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[0]         = big;
  x[(N/2)*incx] = -big;

  res = rzsum(N, x, incx);
  if (res != ref) {
    printf("rzsum(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  res = rdzasum(N, x, incx);
  if (res != refa) {
    printf("rdzasum(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refa), ZIMAG_(refa));
    return 1;
  }

  res = rzdotu(N, x, incx, y, incy);
  if (res != refdu) {
    printf("rzdotu(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refdu), ZIMAG_(refdu));
    return 1;
  }

  res = rzdotc(N, x, incx, y, incy);
  if (res != refdc) {
    printf("rzdotc(x) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refdc), ZIMAG_(refdc));
    return 1;
  }

  //1 Big pos neg at ends
  util_zvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[0]         = big;
  x[(N-1)*incx] = -big;

  res = rzsum(N, x, incx);
  if (res != ref) {
    printf("rzsum(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  res = rdzasum(N, x, incx);
  if (res != refa) {
    printf("rdzasum(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refa), ZIMAG_(refa));
    return 1;
  }

  res = rzdotu(N, x, incx, y, incy);
  if (res != refdu) {
    printf("rzdotu(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refdu), ZIMAG_(refdu));
    return 1;
  }

  res = rzdotc(N, x, incx, y, incy);
  if (res != refdc) {
    printf("rzdotc(x) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refdc), ZIMAG_(refdc));
    return 1;
  }

  //1 Big pos neg at end
  util_zvec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[(N/2)*incx] = big;
  x[(N-1)*incx] = -big;

  res = rzsum(N, x, incx);
  if (res != ref) {
    printf("rzsum(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
    return 1;
  }

  res = rdzasum(N, x, incx);
  if (res != refa) {
    printf("rdzasum(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refa), ZIMAG_(refa));
    return 1;
  }

  res = rzdotu(N, x, incx, y, incy);
  if (res != refdu) {
    printf("rzdotu(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refdu), ZIMAG_(refdu));
    return 1;
  }

  res = rzdotc(N, x, incx, y, incy);
  if (res != refdc) {
    printf("rzdotc(x) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", ZREAL_(res), ZIMAG_(res), ZREAL_(refdc), ZIMAG_(refdc));
    return 1;
  }

  free(x);
  free(y);

  return 0;
}
