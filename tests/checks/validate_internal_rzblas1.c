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

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
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
  double complex *X    = util_zvec_alloc(N, incX);
  double complex *Y    = util_zvec_alloc(N, incY);

  //fill Y with -i where necessary
  util_zvec_fill(N, Y, incY, util_Vec_Constant, -_Complex_I, 1.0);

  //1 Big
  ref   = (N - 1) * small + big;
  refa  = ref;
  refdu  = ref * -_Complex_I;
  refdc  = ref * _Complex_I;

  //1 Big at beginning
  util_zvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = big;

rzsum_sub(N, X, incX, &  res);
  if (res != ref) {
    printf("rzsum(X) = %g + %gi != %g + %gi (1 Big at beginning)\n", creal(res), cimag(res), creal(ref), cimag(ref));
    return 1;
  }

  res = rdzasum(N, X, incX);
  if (res != refa) {
    printf("rdzasum(X) = %g + %gi != %g + %gi (1 Big at beginning)\n", creal(res), cimag(res), creal(refa), cimag(refa));
    return 1;
  }

rzdotu_sub(N, X, incX, Y, incY, &  res);
  if (res != refdu) {
    printf("rzdotu(X) = %g + %gi != %g + %gi (1 Big at beginning)\n", creal(res), cimag(res), creal(refdu), cimag(refdu));
    return 1;
  }

rzdotc_sub(N, X, incX, Y, incY, &  res);
  if (res != refdc) {
    printf("rzdotc(X) = %g + %gi != %g + %gi (1 Big at beginning)\n", creal(res), cimag(res), creal(refdc), cimag(refdc));
    return 1;
  }

  //1 Big in middle
  util_zvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N/2)*incX] = big;

rzsum_sub(N, X, incX, &  res);
  if (res != ref) {
    printf("rzsum(X) = %g + %gi != %g + %gi (1 Big in middle)\n", creal(res), cimag(res), creal(ref), cimag(ref));
    return 1;
  }

  res = rdzasum(N, X, incX);
  if (res != refa) {
    printf("rdzasum(X) = %g + %gi != %g + %gi (1 Big in middle)\n", creal(res), cimag(res), creal(refa), cimag(refa));
    return 1;
  }

rzdotu_sub(N, X, incX, Y, incY, &  res);
  if (res != refdu) {
    printf("rzdotu(X) = %g + %gi != %g + %gi (1 Big in middle)\n", creal(res), cimag(res), creal(refdu), cimag(refdu));
    return 1;
  }

rzdotc_sub(N, X, incX, Y, incY, &  res);
  if (res != refdc) {
    printf("rzdotc(X) = %g + %gi != %g + %gi (1 Big in middle)\n", creal(res), cimag(res), creal(refdc), cimag(refdc));
    return 1;
  }

  //1 Big at end
  util_zvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N-1)*incX] = big;

rzsum_sub(N, X, incX, &  res);
  if (res != ref) {
    printf("rzsum(X) = %g + %gi != %g + %gi (1 Big at end)\n", creal(res), cimag(res), creal(ref), cimag(ref));
    return 1;
  }

  res = rdzasum(N, X, incX);
  if (res != refa) {
    printf("rdzasum(X) = %g + %gi != %g + %gi (1 Big at end)\n", creal(res), cimag(res), creal(refa), cimag(refa));
    return 1;
  }

rzdotu_sub(N, X, incX, Y, incY, &  res);
  if (res != refdu) {
    printf("rzdotu(X) = %g + %gi != %g + %gi (1 Big at end)\n", creal(res), cimag(res), creal(refdu), cimag(refdu));
    return 1;
  }

rzdotc_sub(N, X, incX, Y, incY, &  res);
  if (res != refdc) {
    printf("rzdotc(X) = %g + %gi != %g + %gi (1 Big at end)\n", creal(res), cimag(res), creal(refdc), cimag(refdc));
    return 1;
  }

  //1 Big pos neg
  ref   = (N - 2) * small;
  refa  = ((N - 2) * small) + 2 * big;
  refdu  = ref * -_Complex_I;
  refdc  = ref * _Complex_I;

  //1 Big pos neg at beginning
  util_zvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = big;
  X[(N/2)*incX] = -big;

rzsum_sub(N, X, incX, &  res);
  if (res != ref) {
    printf("rzsum(X) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", creal(res), cimag(res), creal(ref), cimag(ref));
    return 1;
  }

  res = rdzasum(N, X, incX);
  if (res != refa) {
    printf("rdzasum(X) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", creal(res), cimag(res), creal(refa), cimag(refa));
    return 1;
  }

rzdotu_sub(N, X, incX, Y, incY, &  res);
  if (res != refdu) {
    printf("rzdotu(X) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", creal(res), cimag(res), creal(refdu), cimag(refdu));
    return 1;
  }

rzdotc_sub(N, X, incX, Y, incY, &  res);
  if (res != refdc) {
    printf("rzdotc(X) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", creal(res), cimag(res), creal(refdc), cimag(refdc));
    return 1;
  }

  //1 Big pos neg at ends
  util_zvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = big;
  X[(N-1)*incX] = -big;

rzsum_sub(N, X, incX, &  res);
  if (res != ref) {
    printf("rzsum(X) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", creal(res), cimag(res), creal(ref), cimag(ref));
    return 1;
  }

  res = rdzasum(N, X, incX);
  if (res != refa) {
    printf("rdzasum(X) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", creal(res), cimag(res), creal(refa), cimag(refa));
    return 1;
  }

rzdotu_sub(N, X, incX, Y, incY, &  res);
  if (res != refdu) {
    printf("rzdotu(X) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", creal(res), cimag(res), creal(refdu), cimag(refdu));
    return 1;
  }

rzdotc_sub(N, X, incX, Y, incY, &  res);
  if (res != refdc) {
    printf("rzdotc(X) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", creal(res), cimag(res), creal(refdc), cimag(refdc));
    return 1;
  }

  //1 Big pos neg at end
  util_zvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N/2)*incX] = big;
  X[(N-1)*incX] = -big;

rzsum_sub(N, X, incX, &  res);
  if (res != ref) {
    printf("rzsum(X) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", creal(res), cimag(res), creal(ref), cimag(ref));
    return 1;
  }

  res = rdzasum(N, X, incX);
  if (res != refa) {
    printf("rdzasum(X) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", creal(res), cimag(res), creal(refa), cimag(refa));
    return 1;
  }

rzdotu_sub(N, X, incX, Y, incY, &  res);
  if (res != refdu) {
    printf("rzdotu(X) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", creal(res), cimag(res), creal(refdu), cimag(refdu));
    return 1;
  }

rzdotc_sub(N, X, incX, Y, incY, &  res);
  if (res != refdc) {
    printf("rzdotc(X) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", creal(res), cimag(res), creal(refdc), cimag(refdc));
    return 1;
  }

  free(X);
  free(Y);

  return 0;
}
