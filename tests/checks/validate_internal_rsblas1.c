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

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
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
  float *X    = util_svec_alloc(N, incX);
  float *Y    = util_svec_alloc(N, incY);

  //fill Y with -1 where necessary
  util_svec_fill(N, Y, incY, util_Vec_Constant, -1.0, 1.0);

  //1 Big
  ref   = (N - 1) * small + big;
  refa  = ref;
  refd  = ref * -1;

  //1 Big at beginning
  util_svec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = big;

  res = rssum(N, X, incX);
  if (res != ref) {
    printf("rssum(X) = %g != %g (1 Big at beginning)\n", res, ref);
    return 1;
  }

  res = rsasum(N, X, incX);
  if (res != refa) {
    printf("rsasum(X) = %g != %g (1 Big at beginning)\n", res, refa);
    return 1;
  }

  res = rsdot(N, X, incX, Y, incY);
  if (res != refd) {
    printf("rsdot(X) = %g != %g (1 Big at beginning)\n", res, refd);
    return 1;
  }

  //1 Big in middle
  util_svec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N/2)*incX] = big;

  res = rssum(N, X, incX);
  if (res != ref) {
    printf("rssum(X) = %g != %g (1 Big in middle)\n", res, ref);
    return 1;
  }

  res = rsasum(N, X, incX);
  if (res != refa) {
    printf("rsasum(X) = %g != %g (1 Big in middle)\n", res, refa);
    return 1;
  }

  res = rsdot(N, X, incX, Y, incY);
  if (res != refd) {
    printf("rsdot(X) = %g != %g (1 Big in middle)\n", res, refd);
    return 1;
  }

  //1 Big at end
  util_svec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N-1)*incX] = big;

  res = rssum(N, X, incX);
  if (res != ref) {
    printf("rssum(X) = %g != %g (1 Big at end)\n", res, ref);
    return 1;
  }

  res = rsasum(N, X, incX);
  if (res != refa) {
    printf("rsasum(X) = %g != %g (1 Big at end)\n", res, refa);
    return 1;
  }

  res = rsdot(N, X, incX, Y, incY);
  if (res != refd) {
    printf("rsdot(X) = %g != %g (1 Big at end)\n", res, refd);
    return 1;
  }

  //1 Big pos neg
  ref   = (N - 2) * small;
  refa  = ((N - 2) * small) + 2 * big;
  refd  = ref * -1;

  //1 Big pos neg at beginning
  util_svec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = big;
  X[(N/2)*incX] = -big;

  res = rssum(N, X, incX);
  if (res != ref) {
    printf("rssum(X) = %g != %g (1 Big pos neg at beginning)\n", res, ref);
    return 1;
  }

  res = rsasum(N, X, incX);
  if (res != refa) {
    printf("rsasum(X) = %g != %g (1 Big pos neg at beginning)\n", res, refa);
    return 1;
  }

  res = rsdot(N, X, incX, Y, incY);
  if (res != refd) {
    printf("rsdot(X) = %g != %g (1 Big pos neg at beginning)\n", res, refd);
    return 1;
  }

  //1 Big pos neg at ends
  util_svec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = big;
  X[(N-1)*incX] = -big;

  res = rssum(N, X, incX);
  if (res != ref) {
    printf("rssum(X) = %g != %g (1 Big pos neg at ends)\n", res, ref);
    return 1;
  }

  res = rsasum(N, X, incX);
  if (res != refa) {
    printf("rsasum(X) = %g != %g (1 Big pos neg at ends)\n", res, refa);
    return 1;
  }

  res = rsdot(N, X, incX, Y, incY);
  if (res != refd) {
    printf("rsdot(X) = %g != %g (1 Big pos neg at ends)\n", res, refd);
    return 1;
  }

  //1 Big pos neg at end
  util_svec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N/2)*incX] = big;
  X[(N-1)*incX] = -big;

  res = rssum(N, X, incX);
  if (res != ref) {
    printf("rssum(X) = %g != %g (1 Big pos neg at end)\n", res, ref);
    return 1;
  }

  res = rsasum(N, X, incX);
  if (res != refa) {
    printf("rsasum(X) = %g != %g (1 Big pos neg at end)\n", res, refa);
    return 1;
  }

  res = rsdot(N, X, incX, Y, incY);
  if (res != refd) {
    printf("rsdot(X) = %g != %g (1 Big pos neg at end)\n", res, refd);
    return 1;
  }

  free(X);
  free(Y);

  return 0;
}
