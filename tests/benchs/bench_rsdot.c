#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../common/test_vec.h"
#include "../common/test_time.h"
#include "../common/test_perf.h"
#include <rblas.h>
#include <IndexedFP.h>

#include "vecvec_fill_bench_header.h"

#define NAME_SIZE 100
#define FLOP_PER_N 11

extern const char* vecvec_fill_bench_name(int argc, char** argv){
  static char namebuf[NAME_SIZE];
  snprintf(namebuf, NAME_SIZE * sizeof(char), "Benchmark [rsdot]");
  return namebuf;
}

extern int vecvec_fill_bench_test(int argc, char** argv, int N, int incx, int incy, int type, int perf_unit, int trials){
  int rc = 0;
  float res;
  I_float Ires;
  float *x = svec_alloc(N, incx);
  float *y = svec_alloc(N, incy);

  vec_random_seed();

  //fill empty space with random data to check increment
  svec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);
  svec_fill(N * incy, y, 1, vec_fill_RAND, 1.0, 1.0);

  //fill x
  svec_fill(N, x, incx, type, 1.0, opt_read_float(argc, argv, "-c", 1.0));

  //fill y with 1 where necessary
  svec_fill(N, y, incy, vec_fill_CONSTANT, 1.0, 1.0);

  time_tic();
  for(int i = 0; i < trials; i++){
    res = rsdot(N, x, incx, y, incy);
  }
  time_toc();

  printf("%e\n", perf_output(time_read(), N, trials, FLOP_PER_N, perf_unit, perf_prec_SINGLE));

  free(x);
  free(y);
  return rc;
}
