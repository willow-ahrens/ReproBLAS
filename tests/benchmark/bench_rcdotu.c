#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../common/test_vec.h"
#include "../common/test_time.h"
#include <rblas.h>
#include <IndexedFP.h>

#include "vecvec_fill_bench_header.h"

#define NAME_SIZE 100
#define FLOP_PER_N 2

extern const char* vecvec_fill_bench_name(int argc, char** argv){
  static char namebuf[NAME_SIZE];
  snprintf(namebuf, NAME_SIZE * sizeof(char), "Benchmark [rcdot]");
  return namebuf;
}

extern int vecvec_fill_bench_test(int argc, char** argv, int N, int incx, int incy, int type, int unit, int trials){
  int rc = 0;
  float complex res;
  I_float_Complex Ires;
  float complex *x = cvec_alloc(N, incx);
  float complex *y = cvec_alloc(N, incy);

  vec_random_seed();

  //fill empty space with random data to check increment
  cvec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);
  cvec_fill(N * incy, y, 1, vec_fill_RAND, 1.0, 1.0);

  //fill x
  cvec_fill(N, x, incx, type, 1.0, opt_read_float(argc, argv, "-c", 1.0));

  //fill y with -i where necessary
  cvec_fill(N, y, incy, vec_fill_CONSTANT, -_Complex_I, 1.0);

  time_tic();
  for(int i = 0; i < trials; i++){
    res = rcdotu(N, x, incx, y, incy);
  }
  time_toc();

  switch(unit){
    case UNIT_HERTZ:
      printf("%e\n", (double)((long double)N * trials)/time_read());
      break;
    case UNIT_FLOPS:
      printf("%e\n", (double)((long double)N * FLOP_PER_N * trials)/time_read());
      break;
    default:
      rc = 1;
  }

  return rc;
}
