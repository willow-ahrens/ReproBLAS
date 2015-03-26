#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/test_vecvec_header.h"

int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  return "Validate ufpf";
}

int vecvec_test(int argc, char** argv, int N, int incx, int incy) {
  int i;
  float bin = 1.0;
  float ufpx;

  util_random_seed();

  //allocate vector
  float *x = util_svec_alloc(N, incx);

  //check
  for (i = 0; i < N; i++) {
    x[i * incx] = 3 * bin;
    bin *= 2;

    ufpx = ufpf(x[i * incx]);
    if (ufpx != bin) {
      printf("ufpf(%g) = %g != %g\n", x[i * incx], ufpx, bin);
      return 1;
    }
  }

  free(x);

  return 0;
}
