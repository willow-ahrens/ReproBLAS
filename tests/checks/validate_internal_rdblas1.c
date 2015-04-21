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
  return "Validate rdsum, rdasum, rddot (1Big, 1BigPosNeg)";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  double small = 1.0 / (1024.0 * 1024.0);       // 2^-20
  double big   = 1024.0 * 1024.0 * 1024.0 * 32; // 2^35
  double ref;
  double refa;
  double refd;
  double res;

  util_random_seed();

  //allocate vectors
  double *X    = util_dvec_alloc(N, incX);
  double *Y    = util_dvec_alloc(N, incY);

  //fill Y with -1 where necessary
  util_dvec_fill(N, Y, incY, util_Vec_Constant, -1.0, 1.0);

  //1 Big
  ref   = (N - 1) * small + big;
  refa  = ref;
  refd  = ref * -1;

  //1 Big at beginning
  util_dvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = big;

  res = rdsum(N, X, incX);
  if (res != ref) {
    printf("rdsum(X) = %g != %g (1 Big at beginning)\n", res, ref);
    return 1;
  }

  res = rdasum(N, X, incX);
  if (res != refa) {
    printf("rdasum(X) = %g != %g (1 Big at beginning)\n", res, refa);
    return 1;
  }

  res = rddot(N, X, incX, Y, incY);
  if (res != refd) {
    printf("rddot(X) = %g != %g (1 Big at beginning)\n", res, refd);
    return 1;
  }

  //1 Big in middle
  util_dvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N/2)*incX] = big;

  res = rdsum(N, X, incX);
  if (res != ref) {
    printf("rdsum(X) = %g != %g (1 Big in middle)\n", res, ref);
    return 1;
  }

  res = rdasum(N, X, incX);
  if (res != refa) {
    printf("rdasum(X) = %g != %g (1 Big in middle)\n", res, refa);
    return 1;
  }

  res = rddot(N, X, incX, Y, incY);
  if (res != refd) {
    printf("rddot(X) = %g != %g (1 Big in middle)\n", res, refd);
    return 1;
  }

  //1 Big at end
  util_dvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N-1)*incX] = big;

  res = rdsum(N, X, incX);
  if (res != ref) {
    printf("rdsum(X) = %g != %g (1 Big at end)\n", res, ref);
    return 1;
  }

  res = rdasum(N, X, incX);
  if (res != refa) {
    printf("rdasum(X) = %g != %g (1 Big at end)\n", res, refa);
    return 1;
  }

  res = rddot(N, X, incX, Y, incY);
  if (res != refd) {
    printf("rddot(X) = %g != %g (1 Big at end)\n", res, refd);
    return 1;
  }

  //1 Big pos neg
  ref   = (N - 2) * small;
  refa  = ((N - 2) * small) + 2 * big;
  refd  = ref * -1;

  //1 Big pos neg at beginning
  util_dvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = big;
  X[(N/2)*incX] = -big;

  res = rdsum(N, X, incX);
  if (res != ref) {
    printf("rdsum(X) = %g != %g (1 Big pos neg at beginning)\n", res, ref);
    return 1;
  }

  res = rdasum(N, X, incX);
  if (res != refa) {
    printf("rdasum(X) = %g != %g (1 Big pos neg at beginning)\n", res, refa);
    return 1;
  }

  res = rddot(N, X, incX, Y, incY);
  if (res != refd) {
    printf("rddot(X) = %g != %g (1 Big pos neg at beginning)\n", res, refd);
    return 1;
  }

  //1 Big pos neg at ends
  util_dvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = big;
  X[(N-1)*incX] = -big;

  res = rdsum(N, X, incX);
  if (res != ref) {
    printf("rdsum(X) = %g != %g (1 Big pos neg at ends)\n", res, ref);
    return 1;
  }

  res = rdasum(N, X, incX);
  if (res != refa) {
    printf("rdasum(X) = %g != %g (1 Big pos neg at ends)\n", res, refa);
    return 1;
  }

  res = rddot(N, X, incX, Y, incY);
  if (res != refd) {
    printf("rddot(X) = %g != %g (1 Big pos neg at ends)\n", res, refd);
    return 1;
  }

  //1 Big pos neg at end
  util_dvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N/2)*incX] = big;
  X[(N-1)*incX] = -big;

  res = rdsum(N, X, incX);
  if (res != ref) {
    printf("rdsum(X) = %g != %g (1 Big pos neg at end)\n", res, ref);
    return 1;
  }

  res = rdasum(N, X, incX);
  if (res != refa) {
    printf("rdasum(X) = %g != %g (1 Big pos neg at end)\n", res, refa);
    return 1;
  }

  res = rddot(N, X, incX, Y, incY);
  if (res != refd) {
    printf("rddot(X) = %g != %g (1 Big pos neg at end)\n", res, refd);
    return 1;
  }

  free(X);
  free(Y);

  return 0;
}
