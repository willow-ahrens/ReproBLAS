#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include "../common/test_vec.h"

const char* vecvec_name(int argc, char** argv) {
  return "Validate ufpf";
}

int vecvec_check(int argc, char** argv, int N, int incx, int incy) {
  int i;
  float bin = 1.0;
  float ufpx;

  vec_random_seed();

  //allocate vector
  float *x = svec_alloc(N, incx);
  //fill empty space with random data to check increment
  svec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);

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
