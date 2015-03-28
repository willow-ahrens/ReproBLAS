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
  return "Validate rssum, rsasum, rsdot (1Big, 1BigPosNeg)";
}

int vecvec_test(int argc, char** argv, int N, int incx, int incy) {
  (void)argc;
  (void)argv;
  float small = 1.0 / 1024.0; // 2^-10
  float big   = 1024.0 * 32;  // 2^15
  float ref;
  float refa;
  float refd;
  float res;

  util_random_seed();

  //allocate vectors
  float *x    = util_svec_alloc(N, incx);
  float *y    = util_svec_alloc(N, incy);

  //fill y with -1 where necessary
  util_svec_fill(N, y, incy, util_Vec_Constant, -1.0, 1.0);

  //1 Big
  ref   = (N - 1) * small + big;
  refa  = ref;
  refd  = ref * -1;

  //1 Big at beginning
  util_svec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[0]         = big;

  res = rssum(N, x, incx);
  if (res != ref) {
    printf("rssum(x) = %g != %g (1 Big at beginning)\n", res, ref);
    return 1;
  }

  res = rsasum(N, x, incx);
  if (res != refa) {
    printf("rsasum(x) = %g != %g (1 Big at beginning)\n", res, refa);
    return 1;
  }

  res = rsdot(N, x, incx, y, incy);
  if (res != refd) {
    printf("rsdot(x) = %g != %g (1 Big at beginning)\n", res, refd);
    return 1;
  }

  //1 Big in middle
  util_svec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[(N/2)*incx] = big;

  res = rssum(N, x, incx);
  if (res != ref) {
    printf("rssum(x) = %g != %g (1 Big in middle)\n", res, ref);
    return 1;
  }

  res = rsasum(N, x, incx);
  if (res != refa) {
    printf("rsasum(x) = %g != %g (1 Big in middle)\n", res, refa);
    return 1;
  }

  res = rsdot(N, x, incx, y, incy);
  if (res != refd) {
    printf("rsdot(x) = %g != %g (1 Big in middle)\n", res, refd);
    return 1;
  }

  //1 Big at end
  util_svec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
  x[(N-1)*incx] = big;

  res = rssum(N, x, incx);
  if (res != ref) {
    printf("rssum(x) = %g != %g (1 Big at end)\n", res, ref);
    return 1;
  }

  res = rsasum(N, x, incx);
  if (res != refa) {
    printf("rsasum(x) = %g != %g (1 Big at end)\n", res, refa);
    return 1;
  }

  res = rsdot(N, x, incx, y, incy);
  if (res != refd) {
    printf("rsdot(x) = %g != %g (1 Big at end)\n", res, refd);
    return 1;
  }

  //1 Big pos neg
  ref   = (N - 2) * small;
  refa  = ((N - 2) * small) + 2 * big;
  refd  = ref * -1;

  //1 Big pos neg at beginning
  util_svec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
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
  if (res != refd) {
    printf("rsdot(x) = %g != %g (1 Big pos neg at beginning)\n", res, refd);
    return 1;
  }

  //1 Big pos neg at ends
  util_svec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
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
  if (res != refd) {
    printf("rsdot(x) = %g != %g (1 Big pos neg at ends)\n", res, refd);
    return 1;
  }

  //1 Big pos neg at end
  util_svec_fill(N, x, incx, util_Vec_Constant, small, 1.0);
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
  if (res != refd) {
    printf("rsdot(x) = %g != %g (1 Big pos neg at end)\n", res, refd);
    return 1;
  }

  free(x);
  free(y);

  return 0;
}
