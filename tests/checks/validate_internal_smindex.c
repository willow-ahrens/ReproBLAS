#include <binnedBLAS.h>
#include <binned.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../common/test_vecvec_header.h"
#include "../common/test_util.h"


int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  (void)argc;
  (void)argv;
  return "Verify smindex";
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY) {
  (void)argc;
  (void)argv;
  (void)incY;
  int i;
  int index;

  util_random_seed();

  //allocate vector
  float *X = util_svec_alloc(N * ((FLT_MAX_EXP - FLT_MIN_EXP)/SBWIDTH) + 1, incX);

  //check
  for (i = 0; i < N * ((FLT_MAX_EXP - FLT_MIN_EXP)/SBWIDTH) + 1; i++) {
    X[i * incX] = *binned_smbins(i/N) * (1.25/1.5 + 0.5/1.5 * util_drand());
  }
  for (i = 0; i < N * ((FLT_MAX_EXP - FLT_MIN_EXP)/SBWIDTH) + 1; i++) {
    index = binned_smindex(X + i * incX);
    if (*binned_smbins(index)*(1.25/1.5) > X[i * incX]){
      printf("X < 1.25 * 2^(i * W - 1)\n");
      printf("%g < 1.25 * %g\n", X[i * incX], *binned_smbins(index));
      return 1;
    }
    if (*binned_smbins(index)*(1.75/1.5) < X[i * incX]){
      printf("X > 1.75 * 2^(i * W - 1)\n");
      printf("%g > 1.75 * %g\n", X[i * incX], *binned_smbins(index));
      return 1;
    }
  }
  free(X);
  return 0;
}
