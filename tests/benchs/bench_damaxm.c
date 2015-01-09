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

#define FLOP_PER_N 2

extern const char* vecvec_fill_bench_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [damaxm]");
  return name_buffer;
}

extern int vecvec_fill_bench_test(int argc, char** argv, int N, int incx, int incy, int type, int perf_unit, int trials){
  int rc = 0;
  double res;
  I_double Ires;
  double *x = dvec_alloc(N, incx);
  double *y = dvec_alloc(N, incy);

  vec_random_seed();

  //fill empty space with random data to check increment
  dvec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);
  dvec_fill(N * incy, y, 1, vec_fill_RAND, 1.0, 1.0);

  //fill x
  dvec_fill(N, x, incx, type, 1.0, opt_read_float(argc, argv, "-c", 1.0));

  //fill y with 1 where necessary
  dvec_fill(N, y, incy, vec_fill_CONSTANT, 1.0, 1.0);

  time_tic();
  for(int i = 0; i < trials; i++){
    res = damaxm(N, x, incx, y, incy);
  }
  time_toc();

  printf("%e\n", perf_output(time_read(), N, trials, FLOP_PER_N, perf_unit, perf_prec_DOUBLE));

  free(x);
  free(y);
  return rc;
}
