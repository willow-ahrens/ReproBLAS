#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/test_vecvec_header.h"

int vecvec_show_help(void){
  return 0;
}

const char* vecvec_name(int argc, char** argv) {
  return "Validate ufp";
}

int vecvec_test(int argc, char** argv, int N, int incx, int incy) {
  int i;
  double bin = 1.0;
  double ufpx;

  util_random_seed();

  //allocate vector
  double *x = dvec_alloc(N, incx);

  //check
  for (i = 0; i < N; i++) {
    x[i * incx] = 3 * bin;
    bin *= 2;

    ufpx = ufp(x[i * incx]);
    if (ufpx != bin) {
      printf("ufp(%g) = %g != %g\n", x[i * incx], ufpx, bin);
      return 1;
    }
  }
  return 0;
}
