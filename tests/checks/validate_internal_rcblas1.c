#include <idxdBLAS.h>
#include <idxd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../common/test_vecvec_header.h"

int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  (void)argc;
  (void)argv;
  return "Validate rcsum, rscasum, rcdotu/c (+Big, +-Big)";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  float small = 1.0 / (1024.0 * 4.0); // 2^-12
  float big   = 1024.0 * 8.0;         // 2^13
  float complex ref;
  float complex refa;
  float complex refdu;
  float complex refdc;
  float complex res;

  util_random_seed();

  //allocate vectors
  float complex *X    = util_cvec_alloc(N, incX);
  float complex *Y    = util_cvec_alloc(N, incY);

  //fill Y with -i where necessary
  util_cvec_fill(N, Y, incY, util_Vec_Constant, -_Complex_I, 1.0);

  //1 Big
  ref   = (N - 1) * small + big;
  refa  = ref;
  refdu  = ref * -_Complex_I;
  refdc  = ref * _Complex_I;

  //1 Big at beginning
  util_cvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = big;

reproBLAS_rcsum_sub(N, X, incX, &  res);
  if (res != ref) {
    printf("rcsum(X) = %g + %gi != %g + %gi (1 Big at beginning)\n", crealf(res), cimagf(res), crealf(ref), cimagf(ref));
    return 1;
  }

  res = reproBLAS_rscasum(N, X, incX);
  if (res != refa) {
    printf("reproBLAS_rscasum(X) = %g + %gi != %g + %gi (1 Big at beginning)\n", crealf(res), cimagf(res), crealf(refa), cimagf(refa));
    return 1;
  }

reproBLAS_rcdotu_sub(N, X, incX, Y, incY, &  res);
  if (res != refdu) {
    printf("rcdotu(X) = %g + %gi != %g + %gi (1 Big at beginning)\n", crealf(res), cimagf(res), crealf(refdu), cimagf(refdu));
    return 1;
  }

reproBLAS_rcdotc_sub(N, X, incX, Y, incY, &  res);
  if (res != refdc) {
    printf("rcdotc(X) = %g + %gi != %g + %gi (1 Big at beginning)\n", crealf(res), cimagf(res), crealf(refdc), cimagf(refdc));
    return 1;
  }

  //1 Big in middle
  util_cvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N/2)*incX] = big;

reproBLAS_rcsum_sub(N, X, incX, &  res);
  if (res != ref) {
    printf("rcsum(X) = %g + %gi != %g + %gi (1 Big in middle)\n", crealf(res), cimagf(res), crealf(ref), cimagf(ref));
    return 1;
  }

  res = reproBLAS_rscasum(N, X, incX);
  if (res != refa) {
    printf("reproBLAS_rscasum(X) = %g + %gi != %g + %gi (1 Big in middle)\n", crealf(res), cimagf(res), crealf(refa), cimagf(refa));
    return 1;
  }

reproBLAS_rcdotu_sub(N, X, incX, Y, incY, &  res);
  if (res != refdu) {
    printf("rcdotu(X) = %g + %gi != %g + %gi (1 Big in middle)\n", crealf(res), cimagf(res), crealf(refdu), cimagf(refdu));
    return 1;
  }

reproBLAS_rcdotc_sub(N, X, incX, Y, incY, &  res);
  if (res != refdc) {
    printf("rcdotc(X) = %g + %gi != %g + %gi (1 Big in middle)\n", crealf(res), cimagf(res), crealf(refdc), cimagf(refdc));
    return 1;
  }

  //1 Big at end
  util_cvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N-1)*incX] = big;

reproBLAS_rcsum_sub(N, X, incX, &  res);
  if (res != ref) {
    printf("rcsum(X) = %g + %gi != %g + %gi (1 Big at end)\n", crealf(res), cimagf(res), crealf(ref), cimagf(ref));
    return 1;
  }

  res = reproBLAS_rscasum(N, X, incX);
  if (res != refa) {
    printf("reproBLAS_rscasum(X) = %g + %gi != %g + %gi (1 Big at end)\n", crealf(res), cimagf(res), crealf(refa), cimagf(refa));
    return 1;
  }

reproBLAS_rcdotu_sub(N, X, incX, Y, incY, &  res);
  if (res != refdu) {
    printf("rcdotu(X) = %g + %gi != %g + %gi (1 Big at end)\n", crealf(res), cimagf(res), crealf(refdu), cimagf(refdu));
    return 1;
  }

reproBLAS_rcdotc_sub(N, X, incX, Y, incY, &  res);
  if (res != refdc) {
    printf("rcdotc(X) = %g + %gi != %g + %gi (1 Big at end)\n", crealf(res), cimagf(res), crealf(refdc), cimagf(refdc));
    return 1;
  }

  //1 Big pos neg
  ref   = (N - 2) * small;
  refa  = ((N - 2) * small) + 2 * big;
  refdu  = ref * -_Complex_I;
  refdc  = ref * _Complex_I;

  //1 Big pos neg at beginning
  util_cvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = big;
  X[(N/2)*incX] = -big;

reproBLAS_rcsum_sub(N, X, incX, &  res);
  if (res != ref) {
    printf("rcsum(X) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", crealf(res), cimagf(res), crealf(ref), cimagf(ref));
    return 1;
  }

  res = reproBLAS_rscasum(N, X, incX);
  if (res != refa) {
    printf("reproBLAS_rscasum(X) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", crealf(res), cimagf(res), crealf(refa), cimagf(refa));
    return 1;
  }

reproBLAS_rcdotu_sub(N, X, incX, Y, incY, &  res);
  if (res != refdu) {
    printf("rcdotu(X) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", crealf(res), cimagf(res), crealf(refdu), cimagf(refdu));
    return 1;
  }

reproBLAS_rcdotc_sub(N, X, incX, Y, incY, &  res);
  if (res != refdc) {
    printf("rcdotc(X) = %g + %gi != %g + %gi (1 Big pos neg at beginning)\n", crealf(res), cimagf(res), crealf(refdc), cimagf(refdc));
    return 1;
  }

  //1 Big pos neg at ends
  util_cvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[0]         = big;
  X[(N-1)*incX] = -big;

reproBLAS_rcsum_sub(N, X, incX, &  res);
  if (res != ref) {
    printf("rcsum(X) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", crealf(res), cimagf(res), crealf(ref), cimagf(ref));
    return 1;
  }

  res = reproBLAS_rscasum(N, X, incX);
  if (res != refa) {
    printf("reproBLAS_rscasum(X) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", crealf(res), cimagf(res), crealf(refa), cimagf(refa));
    return 1;
  }

reproBLAS_rcdotu_sub(N, X, incX, Y, incY, &  res);
  if (res != refdu) {
    printf("rcdotu(X) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", crealf(res), cimagf(res), crealf(refdu), cimagf(refdu));
    return 1;
  }

reproBLAS_rcdotc_sub(N, X, incX, Y, incY, &  res);
  if (res != refdc) {
    printf("rcdotc(X) = %g + %gi != %g + %gi (1 Big pos neg at ends)\n", crealf(res), cimagf(res), crealf(refdc), cimagf(refdc));
    return 1;
  }

  //1 Big pos neg at end
  util_cvec_fill(N, X, incX, util_Vec_Constant, small, 1.0);
  X[(N/2)*incX] = big;
  X[(N-1)*incX] = -big;

reproBLAS_rcsum_sub(N, X, incX, &  res);
  if (res != ref) {
    printf("rcsum(X) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", crealf(res), cimagf(res), crealf(ref), cimagf(ref));
    return 1;
  }

  res = reproBLAS_rscasum(N, X, incX);
  if (res != refa) {
    printf("reproBLAS_rscasum(X) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", crealf(res), cimagf(res), crealf(refa), cimagf(refa));
    return 1;
  }

reproBLAS_rcdotu_sub(N, X, incX, Y, incY, &  res);
  if (res != refdu) {
    printf("rcdotu(X) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", crealf(res), cimagf(res), crealf(refdu), cimagf(refdu));
    return 1;
  }

reproBLAS_rcdotc_sub(N, X, incX, Y, incY, &  res);
  if (res != refdc) {
    printf("rcdotc(X) = %g + %gi != %g + %gi (1 Big pos neg at end)\n", crealf(res), cimagf(res), crealf(refdc), cimagf(refdc));
    return 1;
  }

  free(X);
  free(Y);

  return 0;
}
