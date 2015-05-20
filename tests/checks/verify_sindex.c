#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/test_vecvec_header.h"
#include "../common/test_util.h"


int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  (void)argc;
  (void)argv;
  return "Verify sindex";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  (void)incY;
  int i;
  int index;

  util_random_seed();

  //allocate vector
  float *X = util_svec_alloc(N*(FLT_MAX_EXP - FLT_MIN_EXP) + 1, incX);

  //check
  for (i = 0; i < (FLT_MAX_EXP - FLT_MIN_EXP) * N; i++) {
    X[i * incX] = ldexpf(0.5 + 0.5 * util_drand48(), (i/N) + FLT_MIN_EXP);
  }
  X[i * incX] = 0.0;
  for (i = 0; i < N * (FLT_MAX_EXP - FLT_MIN_EXP) + 1; i++) {
    index = sindex(X[i * incX]);
    if (sbound(index + 1) / ldexpf(0.75, 24 - siwidth()) > 2 * fabsf(X[i * incX]) * smcompression()){
      printf("2 * |x| !>= 2^(i * W)\n");
      printf("2 * %g !>= %g\n", fabsf(X[i * incX]) * smcompression(), sbound(index + 1)/ldexpf(0.75, 24 - siwidth()));
      return 1;
    }
    if (sbound(index) / ldexpf(0.75, 24 - siwidth()) <= 2 * fabsf(X[i * incX]) * smcompression()){
      printf("2 * |x| !< 2^((i + 1) * W)\n");
      printf("2 * %g !< %g\n", fabsf(X[i * incX]) * smcompression(), sbound(index)/ldexpf(0.75, 24 - siwidth()));
      return 1;
    }
  }
  return 0;
}
